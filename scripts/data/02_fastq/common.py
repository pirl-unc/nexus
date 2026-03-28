import gzip
from random import randint
from typing import List


def create_long_read_fastq_file(
        sequences: List[str],
        output_fastq_file: str,
        num_reads: List[int],
        read_length: int,
        stranded: bool = False
):
    """
    Create a FASTQ (gzipped) file.

    Parameters:
        sequences               :   List of sequences.
        output_fastq_file       :   Output FASTQ file.
        num_reads               :   Number of reads.
        stranded                :   If True, then only the original sequence
                                    is generated for the FASTQ. If False, at 0.5 probability,
                                    the original sequence and the reverse complemented sequence is generated.
    """
    CHUNK_SIZE = 10_000  # flush to gzip every N reads
    quality_cache: dict[int, str] = {}

    with gzip.open(output_fastq_file, 'wt', compresslevel=1) as f:
        buffer = []
        read_idx = 1

        for i, sequence in enumerate(sequences):
            rc_sequence = reverse_complement(sequence)

            for j in range(num_reads[i]):
                sequence_ = sequence if (stranded or j % 2 == 0) else rc_sequence

                if read_length < len(sequence_):
                    start_pos = randint(0, len(sequence_) - read_length)
                    sequence_ = sequence_[start_pos:start_pos + read_length]

                seq_len = len(sequence_)
                if seq_len not in quality_cache:
                    quality_cache[seq_len] = chr(96) * seq_len

                buffer.append(f'@m1000_{read_idx}/ccs\n{sequence_}\n+\n{quality_cache[seq_len]}\n')
                read_idx += 1

                if len(buffer) >= CHUNK_SIZE:
                    f.write(''.join(buffer))
                    buffer.clear()

        if buffer:
            f.write(''.join(buffer))


def create_paired_end_fastq_files(
        sequences: List[str],
        output_r1_fastq_file: str,
        output_r2_fastq_file: str,
        num_reads: List[int],
        read_length: int,
        min_insert_size: int,
        max_insert_size: int
):
    """
    Create FASTQ files.

    Parameters:
        sequences               :   Sequences.
        output_r1_fastq_file    :   Output R1 FASTQ file.
        output_r2_fastq_file    :   Output R2 FASTQ file.
        num_reads               :   Number of reads.
    """
    CHUNK_SIZE = 10_000
    qual = chr(60 + 33) * read_length  # constant — compute once

    with gzip.open(output_r1_fastq_file, 'wt', compresslevel=1) as f1, \
            gzip.open(output_r2_fastq_file, 'wt', compresslevel=1) as f2:

        buf1 = []
        buf2 = []
        read_id = 1

        for i, sequence in enumerate(sequences):
            seq_len = len(sequence)

            # Precompute fallback reads (used when seq too short for any insert)
            fallback_r1 = sequence[:read_length]
            fallback_r2 = reverse_complement(sequence[-read_length:])
            fallback_has_n = 'n' in fallback_r1.lower() or 'n' in fallback_r2.lower()

            for j in range(num_reads[i]):
                while True:
                    insert_size = randint(min_insert_size, max_insert_size)

                    if seq_len < insert_size:
                        if fallback_has_n:
                            continue
                        r1_seq = fallback_r1
                        r2_seq = fallback_r2
                    else:
                        inner_distance = insert_size - 2 * read_length
                        r1_start = randint(0, seq_len - insert_size)
                        r2_start = r1_start + read_length + inner_distance
                        r1_seq = sequence[r1_start:r1_start + read_length]
                        r2_seq = reverse_complement(sequence[r2_start:r2_start + read_length])
                        if 'n' in r1_seq.lower() or 'n' in r2_seq.lower():
                            continue

                    buf1.append(f'@A:1:A:1:1:{read_id}:2 1:N:0:G+A\n{r1_seq}\n+\n{qual}\n')
                    buf2.append(f'@A:1:A:1:1:{read_id}:2 2:N:0:G+A\n{r2_seq}\n+\n{qual}\n')
                    read_id += 1
                    break

                if len(buf1) >= CHUNK_SIZE:
                    f1.write(''.join(buf1))
                    f2.write(''.join(buf2))
                    buf1.clear()
                    buf2.clear()

        if buf1:
            f1.write(''.join(buf1))
            f2.write(''.join(buf2))


def reverse_complement(sequence: str) -> str:
    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'a': 't',
        't': 'a',
        'c': 'g',
        'g': 'c',
        'N': 'N',
        'n': 'n'
    }
    reverse_complement_seq = ''.join(complement[base] for base in reversed(sequence))
    return reverse_complement_seq
