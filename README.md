# **Multi-Omics-Integration**

## **Workflow**

![wm omic-subtypeclassification](https://user-images.githubusercontent.com/42958809/194195264-bbc1271e-0612-4b94-ba19-2f46fececd11.jpg)

![wm omic-1st_filtering](https://user-images.githubusercontent.com/42958809/194182493-ce4330f0-aab7-4fd6-9f52-c5f9481c8be7.png)

![wm omic-2nd_filtering](https://user-images.githubusercontent.com/42958809/194182500-08b84c1b-b713-4acf-a8ef-aaa09f3bddca.png)

## **Environment**

```
<Clone Repo>
cd Multi-omics-intergration
git clone https://github.com/Jin0331/Multi-omics-intergration.git

<GPU>
conda env create --file conda_env_gpu.yaml
conda activate multiomics

<CPU>
conda env create --file conda_env_cpu.yaml
conda activate multiomics-cpu

```

```
<Subgroup Detection>

optional arguments:
  -h, --help            show this help message and exit
  -b BASE, --base BASE  Root Path
  -c CANCER, --cancer CANCER
                        Types of cancer
  -e CYCLE, --cycle CYCLE
                        configuration of the name of output without a filename extension

example : python src/Multi-omics-integration-subgroup.py \
    -b /home/wmbio/WORK/gitworking/Multi-omics-intergration/ \
    -c COAD \
    -e 1000

<Analysis>

optional arguments:
  -h, --help            show this help message and exit
  -b BASE, --base BASE  Root Path
  -c CANCER, --cancer CANCER
                        Types of cancer
  -a CANCER2, --cancer2 CANCER2
                        Types of cancer2
  -d DEA, --dea DEA     DESeq2(deseq2) or EdgeR(edger) or ALL(all)
  -l LOGFC, --logfc LOGFC
                        log2FC Threshold
  -f FDR, --fdr FDR     False dicovery rate Threshold
  -r SEED, --seed SEED  Random Seed

example : python src/Multi-omics-integration-analysis.py \
         -b /home/wmbio/WORK/gitworking/Multi-omics-intergration/ \
         -c COAD \
         -r 331

```

```@ wmbio.co```
