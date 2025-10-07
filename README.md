# MyFirstPipeline
##  The "MyFirstPipeline" Workflow

The overall vision is to have an automated script, using `bash` and/or `awk`, to execute a standard RNAseq analysis workflow from start to finish. The vision is to be as far along as possible towards making it usable with minimal input needed other than executing the program.

### **Step 1: Prepare the Reference Genome**

In order to align any reads, you must first prepare the reference genome for alignment tools. This is a one-time setup, and very important step.

- **Tool**: `bowtie2-build` or `hisat2-build`.
- **Input**: The *Trypanosoma congolense* genome, which is delivered as a collection of FASTA files for different chromosomes or scaffolds.
- **Action**: You then need to combine these individual FASTA files into a single file with the genome. Then, you run the chosen builder tool against this combined file. This results in an assembly of index files (e.g., with the names ending in `.bt2` or `.ht2`).
- **Purpose**: It creates a highly efficient, searchable database of the genome. Searching the unindexed FASTA text file for millions of reads would be very slow; the index makes this alignment step viable.

---

### **Step 2: Quality Control (QC) on Raw Reads**

It is a quality control check on the raw data before you waste your time and computational resources analyzing it. 

- **Tool**: `fastqc`.
- **Input**: The gzipped, paired-end raw sequence data files (`*.fq.gz`).
- **Action**: Your script should loop over all the fastq files provided in the `/fastq` directory and run `fastqc` on each.
- **Purpose**: `fastqc` generates an HTML report that collates a variety of quality measures across your sequencing reads. It checks such things as:
    - **Per Base Sequence Quality**: Are bases properly identified?
- **Adapter Content**: Are there any remaining pieces of sequencing adapters that need to be trimmed away?
    - **Sequence Duplication Levels**: Are there any sequences that are present at elevated levels, maybe due to technical artifacts?
    As you are not required to
*fix* any problems (such as adapter trimming), you have to do the quality check.

---

### **Step 3: Align Reads to the Genome**

This is where you figure out where in the genome each of your short sequencing reads originated from.

- **Tools**: `bowtie2` or `hisat2`, followed by
    `samtools`.
- **Input**: The preprocessed genome index obtained in Step 1 and the paired-end fastq files.
- **Action**:
    1. Align each pair of fastq reads using the aligner (`bowtie2` or `hisat2`). The aligner maps the reads onto the genome index and finds the best hit.
2. The output from the aligner is a **SAM** (Sequence Alignment/Map) file. That is a big text file.
    3. You must pipe this output directly to

        `samtools` and turn the SAM file into a compressed binary format called **BAM** (`.bam`).
4. The BAM file should then be sorted (usually by genomic coordinate) and indexed using
        
        `samtools sort` and `samtools index`, respectively.
- **Purpose**: Alignment gives the raw reads a genomic context. Compressing the alignment data is needed in indexed BAM format for making it compact and allowing for fast data retrieval in the next step.

---

### **Step 4: Quantify Gene Expression (Generate Counts)**

Once you've viewed where the reads overlap, you can count how many reads fall within the boundaries of each gene. This count is applied as a proxy for how active or "expressed" that gene was in the sample.

- **Tool**: `bedtools`.
- **Input**: Each sample's indexed, sorted BAM file from Step 3 and the provided gene annotation file, `TriTrypDB-46_TcongolenseIL3000_2019.bed`.
- **Action**: You feed the alignment file (BAM) and gene location file (BED) into a tool like `bedtools multicov`. The tool will then go through the BED file, gene by gene, and count how many reads in the BAM file overlap with each gene's coordinates.
- **Purpose**: This step transforms the alignment data into a numerical measure of expression for every gene in the genome per sample. The result will be a simple-to-read count table (e.g., Gene ID, Read Count). The gene names can be found in the 4th column of the BED file.

---

### **Step 5 & 6: Aggregate Data and Calculate Fold Change**

These final steps are data manipulation to reach the core biological question: which genes changed their expression levels between different experimental conditions? This is where `awk` really comes into its own. ????

- **Tools**: `awk`, `sort`, `join`, `paste`.
- **Input**: The individual gene count files from Step 4 and the `TriTrypDB-46_TcongolenseIL3000_2019.bed` file to obtain gene descriptions.
- **Actions:**
1. **Aggregate Replicates**: The experiment contains a number of biological replicates across conditions. You will need to utilize
        `awk` to load the `Tco2.fqfiles` sheet  to sum these replicates. Next, calculate
**statistical mean (average)** of the counts for each gene across the replicates in a group (e.g., the average expression of Gene X in the "Clone1 induced 24h" group).
2. **Add Descriptions**: Your output needs to be biologist-friendly. Append or copy the gene descriptions from your BED file onto your table of mean counts.
3. **Calculate Fold Change**: Compare "group-wise" by calculating the fold change (e.g., `mean(group A) / mean(group B)`). This tells you how much a gene's expression has increased or decreased in one condition versus another.
4. **Sort the Output**: The final fold-change data needs to be sorted in **decreasing order of absolute magnitude**. This puts genes with the largest changes at the beginning of the list, which makes them easy to identify. Again, add gene descriptions to this final output file.
