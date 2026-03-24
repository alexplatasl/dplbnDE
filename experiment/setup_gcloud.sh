#!/bin/bash
# =============================================================================
# Google Cloud VM setup and experiment execution
# Creates a 32-core VM, runs the experiment, downloads results, destroys VM
#
# Prerequisites:
#   1. Install Google Cloud SDK: https://cloud.google.com/sdk/docs/install
#   2. Run: gcloud auth login
#   3. Run: gcloud config set project YOUR_PROJECT_ID
#
# Usage: bash experiment/setup_gcloud.sh
# =============================================================================

set -euo pipefail

# Configuration
PROJECT_ID=$(gcloud config get-value project 2>/dev/null)
ZONE="us-central1-a"
VM_NAME="dplbnde-experiment"
MACHINE_TYPE="c2-standard-32"  # 32 vCPUs, 128 GB RAM (~$1.36/hr)
IMAGE_FAMILY="ubuntu-2204-lts"
IMAGE_PROJECT="ubuntu-os-cloud"
DISK_SIZE="50GB"

echo "=== dplbnDE Experiment Runner ==="
echo "Project: $PROJECT_ID"
echo "VM: $VM_NAME ($MACHINE_TYPE)"
echo "Estimated cost: ~\$3-6 USD total"
echo ""

# Step 1: Create VM
echo "[1/7] Creating VM..."
gcloud compute instances create $VM_NAME \
  --zone=$ZONE \
  --machine-type=$MACHINE_TYPE \
  --image-family=$IMAGE_FAMILY \
  --image-project=$IMAGE_PROJECT \
  --boot-disk-size=$DISK_SIZE \
  --boot-disk-type=pd-ssd \
  --scopes=default \
  --quiet

echo "Waiting for VM to be ready..."
sleep 30

# Step 2: Install R and dependencies
echo "[2/7] Installing R..."
gcloud compute ssh $VM_NAME --zone=$ZONE --command="
  sudo apt-get update -qq
  sudo apt-get install -y -qq r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev
  sudo Rscript -e 'install.packages(c(\"devtools\", \"parallel\", \"matrixStats\"), repos=\"https://cran.r-project.org\", quiet=TRUE)'
  sudo Rscript -e 'install.packages(\"bnclassify\", repos=\"https://cran.r-project.org\", quiet=TRUE)'
"

# Step 3: Install dplbnDE from GitHub
echo "[3/7] Installing dplbnDE..."
gcloud compute ssh $VM_NAME --zone=$ZONE --command="
  sudo Rscript -e 'devtools::install_github(\"alexplatasl/dplbnDE\", quiet=TRUE)'
"

# Step 4: Upload data and experiment scripts
echo "[4/7] Uploading data and scripts..."
gcloud compute scp --zone=$ZONE --recurse \
  data/ $VM_NAME:~/data/
gcloud compute scp --zone=$ZONE --recurse \
  experiment/ $VM_NAME:~/experiment/

# Step 5: Run experiment
echo "[5/7] Running experiment (this will take 2-4 hours)..."
gcloud compute ssh $VM_NAME --zone=$ZONE --command="
  cd ~
  export NCORES=32
  export DATA_DIR=~/data
  export OUT_DIR=~/experiment/results
  nohup Rscript experiment/run_experiment.R > experiment/experiment.log 2>&1
  echo 'Experiment finished!'
  tail -5 experiment/experiment.log
"

# Step 6: Run statistical analysis
echo "[6/7] Running statistical analysis..."
gcloud compute ssh $VM_NAME --zone=$ZONE --command="
  cd ~
  export OUT_DIR=~/experiment/results
  Rscript experiment/statistical_tests.R > experiment/stats.log 2>&1
  echo 'Analysis finished!'
  cat experiment/stats.log
"

# Step 7: Download results
echo "[7/7] Downloading results..."
mkdir -p experiment/results
gcloud compute scp --zone=$ZONE --recurse \
  $VM_NAME:~/experiment/results/ experiment/results/
gcloud compute scp --zone=$ZONE \
  $VM_NAME:~/experiment/experiment.log experiment/
gcloud compute scp --zone=$ZONE \
  $VM_NAME:~/experiment/stats.log experiment/

echo ""
echo "=== Results downloaded to experiment/results/ ==="
echo ""

# Destroy VM
read -p "Destroy VM to stop billing? (y/n) " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
  gcloud compute instances delete $VM_NAME --zone=$ZONE --quiet
  echo "VM destroyed."
else
  echo "VM still running at $VM_NAME. Remember to delete it!"
  echo "  gcloud compute instances delete $VM_NAME --zone=$ZONE"
fi
