import os
import glob

# Define the directory containing the FASTQ files
directory = "C:/University/2025/kyle_research_project/Data_processing/Seq_data/master_F_rbcl/trimmed_rbcl"

# Define the combined output manifest file path
combined_manifest_file = os.path.join(directory, "manifest_rbcl.tsv")

# Define sample IDs for each group
impala_samples = [f"{i}-rbcl" for i in range(25, 49)] + ["C1-rbcl"]
elephant_samples = [f"{i}-rbcl" for i in range(49, 73)] + ["C3-rbcl"]

# Convert expected sample IDs to lowercase for case-insensitive comparison
impala_samples_lower = [s.lower() for s in impala_samples]
elephant_samples_lower = [s.lower() for s in elephant_samples]

# Debug: Print the directory being searched
print(f"Searching for .fastq.gz files in: {directory}")

# Get a list of all .fastq.gz files in the directory (exclude .cutadapt.log files)
fastq_files = glob.glob(os.path.join(directory, "*.fastq.gz"))

# Debug: Print the number of files found and their names
print(f"Found {len(fastq_files)} .fastq.gz files:")
for f in fastq_files:
    print(f" - {os.path.basename(f)}")

# Function to convert Windows path to Linux (WSL) path
def to_linux_path(windows_path):
    abs_path = os.path.abspath(windows_path)
    if abs_path.startswith("C:\\"):
        linux_path = abs_path.replace("C:\\", "/mnt/c/").replace("\\", "/")
    else:
        linux_path = abs_path.replace("\\", "/")
    return linux_path

# Function to write a combined manifest (Impala first, then Elephant)
def write_combined_manifest(file_path):
    with open(file_path, "w", newline="\n") as f:
        f.write("sample-id\tabsolute-filepath\n")

        # Process Impala first
        for group_samples, group_samples_lower in [
            (impala_samples, impala_samples_lower),
            (elephant_samples, elephant_samples_lower)
        ]:
            for fq_path in sorted(fastq_files):
                file_name = os.path.basename(fq_path)
                sample_id = file_name.split("_")[0].strip()
                sample_id_lower = sample_id.lower()

                if sample_id_lower in group_samples_lower:
                    original_sample_id = group_samples[group_samples_lower.index(sample_id_lower)]
                    linux_path = to_linux_path(fq_path)
                    f.write(f"{original_sample_id}\t{linux_path}\n")

# Create the combined manifest
write_combined_manifest(combined_manifest_file)
print(f"Combined manifest file created: {combined_manifest_file}")

# Check for missing samples
all_expected_samples = impala_samples + elephant_samples
all_expected_samples_lower = [s.lower() for s in all_expected_samples]
found_samples = [os.path.basename(f).split("_")[0].strip().lower() for f in fastq_files]
missing_samples = [s for s, s_lower in zip(all_expected_samples, all_expected_samples_lower) if s_lower not in found_samples]
if missing_samples:
    print("Warning: The following expected samples were not found in the directory:")
    for sample in missing_samples:
        print(f" - {sample}")
else:
    print("All expected samples were found.")
