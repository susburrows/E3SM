#!/usr/bin/env bash

set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
prepare_script="$script_dir/prepare_locked_dependencies.sh"
fixture_root=$(mktemp -d /private/tmp/ndrop-lock-test.XXXXXX)
trap 'rm -rf -- "$fixture_root"' EXIT

new_repo() {
  local name=$1
  local path="$fixture_root/sources/$name"
  mkdir -p "$path"
  git -C "$path" init -q
  git -C "$path" config user.name "ndrop lock test"
  git -C "$path" config user.email "ndrop-lock-test@example.invalid"
  printf '%s\n' "$name" >"$path/content.txt"
  git -C "$path" add content.txt
  git -C "$path" commit -q -m "initialize $name"
  git -C "$path" rev-parse HEAD
}

mkdir -p "$fixture_root/sources"
ekat_sha=$(new_repo ekat)
kokkos_sha=$(new_repo kokkos)
mam4xx_sha=$(new_repo mam4xx)
e3sm_sha=$(new_repo e3sm)

git -C "$fixture_root/sources/ekat" update-index --add --cacheinfo \
  "160000,$kokkos_sha,extern/kokkos"
git -C "$fixture_root/sources/ekat" commit -q -m "lock kokkos"
ekat_sha=$(git -C "$fixture_root/sources/ekat" rev-parse HEAD)

git -C "$fixture_root/sources/e3sm" update-index --add --cacheinfo \
  "160000,$ekat_sha,externals/ekat"
git -C "$fixture_root/sources/e3sm" update-index --add --cacheinfo \
  "160000,$mam4xx_sha,externals/mam4xx"
git -C "$fixture_root/sources/e3sm" commit -q -m "lock dependencies"
e3sm_sha=$(git -C "$fixture_root/sources/e3sm" rev-parse HEAD)

lock_file="$fixture_root/dependency_lock.json"
python3 - "$lock_file" "$fixture_root/sources" "$e3sm_sha" "$ekat_sha" \
  "$kokkos_sha" "$mam4xx_sha" <<'PY'
import json
import os
import sys

path, sources, e3sm, ekat, kokkos, mam4xx = sys.argv[1:]
repositories = {
    "e3sm": {
        "url": os.path.join(sources, "e3sm"),
        "sha": e3sm,
        "gitlinks": {
            "externals/ekat": ekat,
            "externals/mam4xx": mam4xx,
        },
    },
    "ekat": {
        "url": os.path.join(sources, "ekat"),
        "sha": ekat,
        "gitlinks": {"extern/kokkos": kokkos},
    },
    "kokkos": {"url": os.path.join(sources, "kokkos"), "sha": kokkos},
    "mam4xx": {"url": os.path.join(sources, "mam4xx"), "sha": mam4xx},
}
with open(path, "w", encoding="utf-8") as stream:
    json.dump(
        {"schema_version": 1, "stacks": {"current": {"repositories": repositories}}},
        stream,
        indent=2,
    )
    stream.write("\n")
PY

prepared="$fixture_root/prepared"
"$prepare_script" --lock "$lock_file" --stack current \
  --source-root "$fixture_root/sources/e3sm" --root "$prepared" \
  --record "$prepared/provenance.json" >/dev/null

test -s "$prepared/provenance.json"
"$prepare_script" --lock "$lock_file" --stack current \
  --source-root "$fixture_root/sources/e3sm" --root "$prepared" \
  --record "$prepared/provenance-verified.json" --verify-only >/dev/null

printf '%s\n' dirty >"$prepared/mam4xx/untracked.txt"
if "$prepare_script" --lock "$lock_file" --stack current \
  --source-root "$fixture_root/sources/e3sm" --root "$prepared" \
  --verify-only >"$fixture_root/dirty.out" 2>&1; then
  echo "dirty checkout was incorrectly accepted" >&2
  exit 1
fi
grep -q "mam4xx checkout is dirty" "$fixture_root/dirty.out"
rm "$prepared/mam4xx/untracked.txt"

python3 - "$lock_file" <<'PY'
import json
import sys

path = sys.argv[1]
with open(path, encoding="utf-8") as stream:
    lock = json.load(stream)
lock["stacks"]["current"]["repositories"]["mam4xx"]["sha"] = "0" * 40
with open(path, "w", encoding="utf-8") as stream:
    json.dump(lock, stream, indent=2)
    stream.write("\n")
PY
if "$prepare_script" --lock "$lock_file" --stack current \
  --source-root "$fixture_root/sources/e3sm" --root "$prepared" \
  --verify-only >"$fixture_root/sha.out" 2>&1; then
  echo "SHA mismatch was incorrectly accepted" >&2
  exit 1
fi
grep -q "mam4xx SHA mismatch" "$fixture_root/sha.out"

echo "RESULT passed=3 failed=0"
