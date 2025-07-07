import os
import sys
import shutil
import subprocess
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed
import pynvml

MINIMAP2 = "nanopore_gpuenv/bin/minimap2"
SAMTOOLS = "nanopore_gpuenv/bin/samtools"
MEDAKA_MODEL = "r1041_e82_400bps_hac_v5.0.0"
DEBUG = False

def get_available_gpus():
    try:
        pynvml.nvmlInit()
        gpu_count = pynvml.nvmlDeviceGetCount()
        return list(range(gpu_count))
    except Exception as e:
        print(f"⚠️ Could not detect GPUs: {e}")
        return [0]

def get_optimal_threads_per_gpu(num_gpus):
    total_cpus = os.cpu_count() or 4
    return max(total_cpus // max(num_gpus, 1), 1)

def process_single_file(args):
    input_path, reference_path, output_dir, min_seqs, gpu_id, threads = args
    fasta_file = os.path.basename(input_path)

    tmp_dir = None  # Initialize to avoid UnboundLocalError

    try:
        seq_count = sum(1 for _ in SeqIO.parse(input_path, "fasta"))
        if seq_count < min_seqs:
            return f"⏭️ Skipped {fasta_file} ({seq_count} sequences)"

        sample_name = os.path.splitext(fasta_file)[0]
        tmp_dir = os.path.join(output_dir, f"tmp_{sample_name}")
        os.makedirs(tmp_dir, exist_ok=True)

        # Alignment steps
        sam_path = os.path.join(tmp_dir, "aligned.sam")
        bam_path = os.path.join(tmp_dir, "aligned.bam")
        sorted_bam = os.path.join(tmp_dir, "aligned.sorted.bam")

        cmd_minimap = [MINIMAP2, "-a", reference_path, input_path]
        minimap_proc = subprocess.run(cmd_minimap, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if minimap_proc.returncode != 0 or not minimap_proc.stdout:
            return f"❌ Minimap2 failed for {fasta_file}: {minimap_proc.stderr.decode()}"

        with open(sam_path, "w") as sam_out:
            sam_out.write(minimap_proc.stdout.decode())

        subprocess.run([SAMTOOLS, "view", "-bS", sam_path, "-o", bam_path], check=True)
        subprocess.run([SAMTOOLS, "sort", bam_path, "-o", sorted_bam], check=True)
        subprocess.run([SAMTOOLS, "index", sorted_bam], check=True)

        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

        cmd_medaka = [
            "nanopore_gpuenv/bin/medaka_consensus",
            "-i", input_path,
            "-d", reference_path,
            "-o", tmp_dir,
            "-t", str(threads),
            "-m", MEDAKA_MODEL
        ]

        medaka_proc = subprocess.run(cmd_medaka, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
        if DEBUG:
            print(f"[{fasta_file}] STDOUT:\n{medaka_proc.stdout}")
        if medaka_proc.returncode != 0:
            return f"❌ Medaka failed for {fasta_file}: {medaka_proc.stderr}"

        consensus_tmp = os.path.join(tmp_dir, "consensus.fasta")
        consensus_out = os.path.join(output_dir, f"{sample_name}_consensus.fa")

        if os.path.exists(consensus_tmp):
            shutil.copyfile(consensus_tmp, consensus_out)
            return f"✅ Done: {fasta_file} → {consensus_out}"
        else:
            return f"❌ consensus.fasta not found for {fasta_file}"

    except Exception as e:
        return f"⚠️ Exception in {fasta_file}: {e}"

    finally:
        if tmp_dir and os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir, ignore_errors=True)

def parallel_consensus(fasta_dir, reference_path, output_dir, min_seqs=2):
    os.makedirs(output_dir, exist_ok=True)
    fasta_files = [
        os.path.join(fasta_dir, f)
        for f in os.listdir(fasta_dir)
        if f.endswith(".fa") or f.endswith(".fasta")
    ]
    print(f"📂 Found {len(fasta_files)} FASTA files in {fasta_dir}")

    gpus = get_available_gpus()
    print(f"🖥️ Detected {len(gpus)} GPU(s): {gpus}")
    threads_per_gpu = get_optimal_threads_per_gpu(len(gpus))
    print(f"🧵 Assigning {threads_per_gpu} thread(s) per GPU process")

    tasks = [
        (fasta_path, reference_path, output_dir, min_seqs, gpus[i % len(gpus)], threads_per_gpu)
        for i, fasta_path in enumerate(fasta_files)
    ]

    with ProcessPoolExecutor(max_workers=len(gpus)) as executor:
        futures = {executor.submit(process_single_file, t): t[0] for t in tasks}
        for future in as_completed(futures):
            print(future.result())

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python generate_consensus_medaka.py <fasta_dir> <reference.fa> <output_dir>")
        sys.exit(1)

    fasta_dir = sys.argv[1]
    reference_path = sys.argv[2]
    output_dir = sys.argv[3]

    parallel_consensus(fasta_dir, reference_path, output_dir)
