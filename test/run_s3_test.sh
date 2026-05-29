#!/usr/bin/env bash
# Run the bgen S3 integration test against a local MinIO container.
#
# Requirements: docker, aws CLI (v2).
# Usage: ./test/run_s3_test.sh [--keep]   (--keep: leave container running after success)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BGEN_FILE="$REPO_ROOT/test/data/example.14bits.bgen"
MINIO_CONTAINER="bgen-minio-test"
MINIO_PORT=19000
BUCKET="bgen-test"
KEY="example.14bits.bgen"
ACCESS="minio_access"
SECRET="minio_secret"
KEEP=0

for arg in "$@"; do
  [[ "$arg" == "--keep" ]] && KEEP=1
done

cleanup() {
  echo "--- Stopping MinIO container ---"
  docker rm -f "$MINIO_CONTAINER" 2>/dev/null || true
}
[[ $KEEP -eq 0 ]] && trap cleanup EXIT

# ── 1. Start MinIO ───────────────────────────────────────────────────────────
echo "--- Starting MinIO ---"
docker rm -f "$MINIO_CONTAINER" 2>/dev/null || true
docker run -d \
  --name "$MINIO_CONTAINER" \
  -p "${MINIO_PORT}:9000" \
  -e "MINIO_ROOT_USER=${ACCESS}" \
  -e "MINIO_ROOT_PASSWORD=${SECRET}" \
  minio/minio server /data

# ── 2. Wait for MinIO to become healthy ──────────────────────────────────────
echo "--- Waiting for MinIO to be ready ---"
for i in $(seq 1 30); do
  if curl -sf "http://127.0.0.1:${MINIO_PORT}/minio/health/live" >/dev/null 2>&1; then
    echo "MinIO is ready."
    break
  fi
  sleep 1
  if [[ $i -eq 30 ]]; then
    echo "MinIO did not become ready in time." >&2
    exit 1
  fi
done

# ── 3. Upload the test BGEN file via AWS CLI ─────────────────────────────────
echo "--- Creating bucket and uploading test file ---"
export AWS_ACCESS_KEY_ID="$ACCESS"
export AWS_SECRET_ACCESS_KEY="$SECRET"
export AWS_DEFAULT_REGION="us-east-1"
export AWS_ENDPOINT_URL="http://127.0.0.1:${MINIO_PORT}"

aws s3 mb "s3://${BUCKET}" --endpoint-url "$AWS_ENDPOINT_URL"
aws s3 cp "$BGEN_FILE" "s3://${BUCKET}/${KEY}" --endpoint-url "$AWS_ENDPOINT_URL"

# Also create a public (anonymous-read) bucket for the no-sign test
BUCKET_PUB="bgen-public"
aws s3 mb "s3://${BUCKET_PUB}" --endpoint-url "$AWS_ENDPOINT_URL"
aws s3api put-bucket-policy \
  --bucket "$BUCKET_PUB" \
  --policy '{"Version":"2012-10-17","Statement":[{"Effect":"Allow","Principal":"*","Action":"s3:GetObject","Resource":"arn:aws:s3:::'"$BUCKET_PUB"'/*"}]}' \
  --endpoint-url "$AWS_ENDPOINT_URL"
aws s3 cp "$BGEN_FILE" "s3://${BUCKET_PUB}/${KEY}" --endpoint-url "$AWS_ENDPOINT_URL"

echo "Uploaded: s3://${BUCKET}/${KEY} (private)"
echo "Uploaded: s3://${BUCKET_PUB}/${KEY} (public)"

# ── 4. Build with S3 support ─────────────────────────────────────────────────
BUILD_DIR="$REPO_ROOT/build-s3-test"
echo "--- Configuring (BGEN_ENABLE_S3=ON) ---"
cmake -B "$BUILD_DIR" "$REPO_ROOT" \
  -DBGEN_ENABLE_S3=ON \
  -DBGEN_BUILD_TESTS=ON \
  -DBGEN_USE_SYSTEM_ZLIB_NG=OFF \
  -DBGEN_USE_SYSTEM_ZSTD=OFF \
  -DBGEN_USE_SYSTEM_ATHR=OFF \
  -DCMAKE_BUILD_TYPE=Debug \
  -DBGEN_ENABLE_SANITIZERS=OFF

echo "--- Building ---"
cmake --build "$BUILD_DIR" --target test_s3_open -- -j"$(sysctl -n hw.logicalcpu 2>/dev/null || nproc)"

# ── 5. Run the tests ──────────────────────────────────────────────────────────
echo "--- Test 1: signed access (explicit credentials) ---"
export BGEN_TEST_S3_URI="s3://${BUCKET}/${KEY}"
export BGEN_TEST_S3_METAFILE="$REPO_ROOT/test/data/example.14bits.bgen.metafile"
export AWS_S3_ADDRESSING_STYLE="path"   # MinIO uses path-style

"$BUILD_DIR/test/test_s3_open"
EXIT_CODE=$?
[[ $EXIT_CODE -ne 0 ]] && { echo "FAILURE: signed test exited $EXIT_CODE" >&2; exit $EXIT_CODE; }

echo ""
echo "--- Test 2: no-sign access (public bucket) ---"
export BGEN_TEST_S3_URI="s3://${BUCKET_PUB}/${KEY}"
unset AWS_ACCESS_KEY_ID
unset AWS_SECRET_ACCESS_KEY
export AWS_NO_SIGN_REQUEST=1

"$BUILD_DIR/test/test_s3_open"
EXIT_CODE=$?
[[ $EXIT_CODE -ne 0 ]] && { echo "FAILURE: no-sign test exited $EXIT_CODE" >&2; exit $EXIT_CODE; }

# ── 6. Python / UPath integration tests ──────────────────────────────────────
echo ""
echo "--- Installing Python package with S3 support ---"

# Install extra test deps into the pixi environment
pixi run -e test pip install "universal_pathlib>=0.2" "s3fs>=2024.1" --quiet

echo "--- Building Python package with S3 ---"
CMAKE_ARGS="-DBGEN_ENABLE_S3=ON" \
  pixi run -e test pip install --no-build-isolation -e "$REPO_ROOT" --quiet

echo ""
echo "--- Test 3: Python UPath credentials ---"
# Unset ambient AWS credentials so tests prove credentials flow from UPath
unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_SESSION_TOKEN \
      AWS_NO_SIGN_REQUEST AWS_ENDPOINT_URL AWS_DEFAULT_REGION

BGEN_TEST_MINIO_ENDPOINT="http://127.0.0.1:${MINIO_PORT}" \
BGEN_TEST_MINIO_ACCESS="$ACCESS" \
BGEN_TEST_MINIO_SECRET="$SECRET" \
BGEN_TEST_S3_METAFILE="$REPO_ROOT/test/data/example.14bits.bgen.metafile" \
  pixi run -e test pytest "$REPO_ROOT/python/tests/test_s3_upath.py" -v
EXIT_CODE=$?

if [[ $EXIT_CODE -eq 0 ]]; then
  echo ""
  echo "SUCCESS: all S3 tests passed."
else
  echo ""
  echo "FAILURE: Python UPath test exited with code $EXIT_CODE." >&2
fi

exit $EXIT_CODE
