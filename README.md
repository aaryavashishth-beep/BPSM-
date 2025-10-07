# MyFirstPipeline 
##  The "MyFirstPipeline" Workflow

The core task is to build an automated script, using `bash` and/or `awk`, that executes a standard RNAseq analysis workflow from start to finish. The goal is to make it user-friendly, requiring minimal input other than executing the program.

### **Step 1: Prepare the Reference Genome**

Before you can align any reads, you must prepare the reference genome for the alignment software. This is a crucial, one-time setup step.

- **Tool**: `bowtie2-build` or `hisat2-build`.
- **Input**: The *Trypanosoma congolense* genome, which is provided as a collection of FASTA files representing different chromosomes or scaffolds.
- **Action**: You need to combine these individual FASTA files into a single genome file. Then, you run the chosen builder tool on this combined file. This creates a set of index files (e.g., with extensions like `.bt2` or `.ht2`).
- **Purpose**: This indexing process creates a highly efficient, searchable database of the genome. Trying to search the raw FASTA text file for millions of reads would be incredibly slow; the index makes the alignment step feasible.

---

### **Step 2: Quality Control (QC) on Raw Reads**

This step is about checking the quality of the raw data before you invest time and computational resources in analyzing it. 🧐

- **Tool**: `fastqc`.
- **Input**: The compressed, paired-end raw sequence data files (`*.fq.gz`).
- **Action**: Your script should loop through all the fastq files provided in the `/fastq` directory and run `fastqc` on each one.
- **Purpose**: `fastqc` generates an HTML report that summarizes various quality metrics for your sequencing reads. It checks things like:
    - **Per Base Sequence Quality**: Are the bases accurately identified?
    - **Adapter Content**: Are there leftover bits of sequencing adapters that need to be trimmed?
    - **Sequence Duplication Levels**: Are some sequences overrepresented, possibly due to technical artifacts?
    While you are not required to
        
        *fix* any issues (like trimming adapters), you must perform the quality check.
        

---

### **Step 3: Align Reads to the Genome**

This is where you determine where in the genome each of your short sequencing reads came from.

- **Tools**: `bowtie2` or `hisat2` , followed by
    
    `samtools`.
    
- **Input**: The prepared genome index from Step 1 and the paired-end fastq files.
- **Action**:
    1. Run the aligner (`bowtie2` or `hisat2`) for each pair of fastq files. The aligner takes the reads and compares them against the genome index to find the best match.
    2. The output of the aligner is a **SAM** (Sequence Alignment/Map) file. This is a very large text file.
    3. You must immediately pipe this output to
        
        `samtools` to convert the SAM file into a compressed binary format called **BAM** (`.bam`).
        
    4. The BAM file should then be sorted (usually by genomic coordinate) and indexed using
        
        `samtools sort` and `samtools index`, respectively.
        
- **Purpose**: Alignment gives the raw reads a genomic context. Converting to an indexed BAM format is essential for making the alignment data compact and allowing for rapid data retrieval in the next step.

---

### **Step 4: Quantify Gene Expression (Generate Counts)**

Now that you know where the reads map, you can count how many reads fall within the boundaries of each gene. This count is used as a proxy for how active or "expressed" that gene was in the sample.

- **Tool**: `bedtools`.
- **Input**: The indexed, sorted BAM file for each sample from Step 3 and the provided gene annotation file, `TriTrypDB-46_TcongolenseIL3000_2019.bed`.
- **Action**: Using a tool like `bedtools multicov`, you provide the alignment file (BAM) and the gene location file (BED). The program will go through the BED file, one gene at a time, and count how many reads in the BAM file overlap with that gene's coordinates.
- **Purpose**: This step translates the alignment information into a quantitative measurement of expression for every gene in the genome for each sample. The output will be a simple count table (e.g., Gene ID, Read Count). The gene names are in the 4th column of the BED file.

---

### **Step 5 & 6: Aggregate Data and Calculate Fold Change**

These final steps involve data manipulation to get to the key biological question: which genes changed their expression levels between different experimental conditions? This is where `awk` becomes particularly powerful. 📊

- **Tools**: `awk`, `sort`, `join`, `paste`.
- **Input**: The individual gene count files from Step 4 and the `TriTrypDB-46_TcongolenseIL3000_2019.bed` file for gene descriptions.
- **Actions**:
    1. **Aggregate Replicates**: The experiment has multiple biological replicates for each condition. You'll need to use
        
        `awk` to parse the `Tco2.fqfiles` sheet  to group these replicates. Then, calculate the
        
        **statistical mean (average)** of the counts for each gene across the replicates within a group (e.g., the average expression of Gene X in the "Clone1 induced 24h" group).
        
    2. **Add Descriptions**: Your output must be biologist-friendly. Join or paste the gene descriptions from the BED file with your table of mean counts.
    3. **Calculate Fold Change**: Perform "group-wise" comparisons by calculating the fold change (e.g., `mean(group A) / mean(group B)`). This tells you how much a gene's expression has increased or decreased in one condition relative to another.
    4. **Sort the Output**: The final fold-change data must be sorted in **decreasing order of absolute magnitude**. This brings the genes with the most dramatic changes to the top of the list, making them easy to identify. Again, include gene descriptions in this final output file.
