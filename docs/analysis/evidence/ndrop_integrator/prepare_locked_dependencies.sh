#!/usr/bin/env bash

set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
lock_file="$script_dir/dependency_lock.json"
stack=current
prepare_root=
results_root="/Users/burr114/e3sm-results/ndrop-screening"
source_root=$(git -C "$script_dir" rev-parse --show-toplevel)
record_file=
verify_only=false
cmake_build_type=Release
cmake_generator="Unix Makefiles"
cmake_flags=
invocation=("$0" "$@")

usage() {
  cat <<'EOF'
Usage: prepare_locked_dependencies.sh [options]

Options:
  --stack NAME          Lock stack to prepare (default: current).
  --root PATH           Destination below /Users/burr114/e3sm-scratch.
  --results-root PATH   Retained output root below /Users/burr114/e3sm-results.
  --source-root PATH    E3SM repository containing the locked base commit.
  --lock PATH           Alternate lock manifest (used by behavior tests).
  --record PATH         Provenance JSON output below the retained-results root.
  --verify-only         Do not clone or fetch; only verify existing checkouts.
  --build-type TYPE     CMake build type recorded in provenance.
  --generator NAME      CMake generator recorded in provenance.
  --cmake-flags TEXT    Exact additional CMake flags recorded in provenance.
  -h, --help            Show this help.

Existing mismatched, attached, or dirty checkouts are never reset or cleaned;
the script fails closed and leaves them for inspection.
EOF
}

while (($#)); do
  case "$1" in
  --stack=*)
    stack=${1#*=}
    shift
    ;;
  --stack)
    stack=$2
    shift 2
    ;;
  --root=*)
    prepare_root=${1#*=}
    shift
    ;;
  --root)
    prepare_root=$2
    shift 2
    ;;
  --results-root=*)
    results_root=${1#*=}
    shift
    ;;
  --results-root)
    results_root=$2
    shift 2
    ;;
  --source-root=*)
    source_root=${1#*=}
    shift
    ;;
  --source-root)
    source_root=$2
    shift 2
    ;;
  --lock=*)
    lock_file=${1#*=}
    shift
    ;;
  --lock)
    lock_file=$2
    shift 2
    ;;
  --record=*)
    record_file=${1#*=}
    shift
    ;;
  --record)
    record_file=$2
    shift 2
    ;;
  --verify-only)
    verify_only=true
    shift
    ;;
  --build-type=*)
    cmake_build_type=${1#*=}
    shift
    ;;
  --build-type)
    cmake_build_type=$2
    shift 2
    ;;
  --generator=*)
    cmake_generator=${1#*=}
    shift
    ;;
  --generator)
    cmake_generator=$2
    shift 2
    ;;
  --cmake-flags=*)
    cmake_flags=${1#*=}
    shift
    ;;
  --cmake-flags)
    cmake_flags=$2
    shift 2
    ;;
  -h | --help)
    usage
    exit 0
    ;;
  *)
    echo "unknown argument: $1" >&2
    usage >&2
    exit 2
    ;;
  esac
done

prepare_root=${prepare_root:-"/Users/burr114/e3sm-scratch/ndrop-screening/dependencies/$stack"}

case "$prepare_root" in
/Users/burr114/e3sm-scratch/*) ;;
*)
  echo "preparation root must be a specific path below /Users/burr114/e3sm-scratch" >&2
  exit 2
  ;;
esac
case "$results_root" in
/Users/burr114/e3sm-results/*) ;;
*)
  echo "results root must be a specific path below /Users/burr114/e3sm-results" >&2
  exit 2
  ;;
esac

if [[ ! -f "$lock_file" ]]; then
  echo "lock manifest not found: $lock_file" >&2
  exit 2
fi
if ! git -C "$source_root" rev-parse --show-toplevel >/dev/null 2>&1; then
  echo "source root is not a Git repository: $source_root" >&2
  exit 2
fi

record_file=${record_file:-"$results_root/provenance/dependencies-$stack.json"}
case "$record_file" in
/Users/burr114/e3sm-results/*) ;;
*)
  echo "provenance record must be below /Users/burr114/e3sm-results" >&2
  exit 2
  ;;
esac
mkdir -p "$prepare_root"

repository_rows=()
while IFS= read -r row; do
  repository_rows+=("$row")
done < <(
  python3 - "$lock_file" "$stack" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as stream:
    lock = json.load(stream)
try:
    repositories = lock["stacks"][sys.argv[2]]["repositories"]
except KeyError as error:
    raise SystemExit(f"missing lock entry: {error}")
for name in sorted(repositories):
    repo = repositories[name]
    print(f"{name}\t{repo['url']}\t{repo['sha']}")
PY
)

if ((${#repository_rows[@]} == 0)); then
  echo "lock stack contains no repositories: $stack" >&2
  exit 2
fi

locked_sha() {
  local requested=$1
  local row name url sha
  for row in "${repository_rows[@]}"; do
    IFS=$'\t' read -r name url sha <<<"$row"
    if [[ "$name" == "$requested" ]]; then
      printf '%s\n' "$sha"
      return 0
    fi
  done
  return 1
}

verify_checkout() {
  local name=$1
  local expected=$2
  local target="$prepare_root/$name"
  local actual branch dirty

  if ! git -C "$target" rev-parse --git-dir >/dev/null 2>&1; then
    echo "$name checkout is not a Git repository: $target" >&2
    return 1
  fi
  actual=$(git -C "$target" rev-parse HEAD)
  if [[ "$actual" != "$expected" ]]; then
    echo "$name SHA mismatch: expected $expected, found $actual" >&2
    return 1
  fi
  branch=$(git -C "$target" symbolic-ref -q --short HEAD || true)
  if [[ -n "$branch" ]]; then
    echo "$name checkout is attached to branch $branch; detached HEAD required" >&2
    return 1
  fi
  dirty=$(git -C "$target" status --porcelain=v1 --untracked-files=all)
  if [[ -n "$dirty" ]]; then
    echo "$name checkout is dirty: $target" >&2
    return 1
  fi
}

prepare_checkout() {
  local name=$1
  local url=$2
  local sha=$3
  local target="$prepare_root/$name"
  local clone_url=$url

  if [[ "$name" == e3sm ]]; then
    clone_url=$source_root
    if ! git -C "$source_root" cat-file -e "$sha^{commit}" 2>/dev/null; then
      echo "locked E3SM commit is absent from source repository: $sha" >&2
      return 1
    fi
  fi

  if [[ ! -e "$target" ]]; then
    if $verify_only; then
      echo "$name checkout is missing in verify-only mode: $target" >&2
      return 1
    fi
    if [[ "$name" == e3sm ]]; then
      # The E3SM checkout is a temporary immutable base used to verify
      # gitlinks. Share the local object store instead of duplicating it.
      git clone --shared --no-checkout "$clone_url" "$target"
    else
      git clone --no-checkout --filter=blob:none "$clone_url" "$target"
    fi
    if ! git -C "$target" cat-file -e "$sha^{commit}" 2>/dev/null; then
      git -C "$target" fetch --no-tags --depth=1 origin "$sha"
    fi
    git -C "$target" checkout --detach "$sha"
  fi
  verify_checkout "$name" "$sha"
}

for row in "${repository_rows[@]}"; do
  IFS=$'\t' read -r name url sha <<<"$row"
  prepare_checkout "$name" "$url" "$sha"
done

gitlink_rows=()
while IFS= read -r row; do
  gitlink_rows+=("$row")
done < <(
  python3 - "$lock_file" "$stack" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as stream:
    repositories = json.load(stream)["stacks"][sys.argv[2]]["repositories"]
for name in sorted(repositories):
    for path, sha in sorted(repositories[name].get("gitlinks", {}).items()):
        print(f"{name}\t{path}\t{sha}")
PY
)

for row in "${gitlink_rows[@]}"; do
  IFS=$'\t' read -r parent path expected <<<"$row"
  actual=$(git -C "$prepare_root/$parent" ls-tree HEAD -- "$path" |
    awk '$1 == "160000" {print $3}')
  if [[ "$actual" != "$expected" ]]; then
    echo "$parent gitlink mismatch at $path: expected $expected, found ${actual:-missing}" >&2
    exit 1
  fi
  child_name=$(basename "$path")
  if child_sha=$(locked_sha "$child_name" 2>/dev/null); then
    if [[ "$child_sha" != "$expected" ]]; then
      echo "lock is internally inconsistent for $parent:$path" >&2
      exit 1
    fi
  fi
done

mkdir -p "$(dirname "$record_file")"
python3 - "$lock_file" "$stack" "$prepare_root" "$source_root" \
  "$record_file" "$cmake_build_type" "$cmake_generator" "$cmake_flags" \
  "${invocation[@]}" <<'PY'
import json
import os
import platform
import shlex
import subprocess
import sys


def output(*command):
    return subprocess.check_output(command, text=True, stderr=subprocess.STDOUT).strip()


def optional_output(*command):
    try:
        return output(*command)
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


(
    lock_path,
    stack,
    prepare_root,
    source_root,
    record_path,
    build_type,
    generator,
    cmake_flags,
    *invocation,
) = sys.argv[1:]
with open(lock_path, encoding="utf-8") as stream:
    lock = json.load(stream)

repositories = {}
for name, entry in lock["stacks"][stack]["repositories"].items():
    path = os.path.join(prepare_root, name)
    status = output("git", "-C", path, "status", "--porcelain=v1", "--untracked-files=all")
    repositories[name] = {
        "path": path,
        "url": entry["url"],
        "locked_sha": entry["sha"],
        "actual_sha": output("git", "-C", path, "rev-parse", "HEAD"),
        "detached": subprocess.run(
            ["git", "-C", path, "symbolic-ref", "-q", "HEAD"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        ).returncode != 0,
        "dirty": bool(status),
    }

source_status = output(
    "git", "-C", source_root, "status", "--porcelain=v1", "--untracked-files=all"
)
cxx = os.environ.get("CXX", "c++")
record = {
    "schema_version": 1,
    "stack": stack,
    "lock_manifest": os.path.abspath(lock_path),
    "invocation": shlex.join(invocation),
    "source_worktree": {
        "path": source_root,
        "head": output("git", "-C", source_root, "rev-parse", "HEAD"),
        "branch": optional_output("git", "-C", source_root, "branch", "--show-current"),
        "dirty": bool(source_status),
    },
    "repositories": repositories,
    "build": {
        "cmake_generator": generator,
        "cmake_build_type": build_type,
        "cmake_flags": cmake_flags,
        "cmake_version": optional_output("cmake", "--version").splitlines()[0],
        "compiler_command": cxx,
        "compiler_version": optional_output(cxx, "--version").splitlines()[0],
    },
    "machine": {
        "os": platform.platform(),
        "architecture": platform.machine(),
        "processor": platform.processor() or "unavailable",
        "cpu": optional_output("sysctl", "-n", "machdep.cpu.brand_string"),
    },
}
with open(record_path, "w", encoding="utf-8") as stream:
    json.dump(record, stream, indent=2, sort_keys=True)
    stream.write("\n")
print(record_path)
PY

echo "locked dependency stack verified: $stack"
echo "provenance record: $record_file"
