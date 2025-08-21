# ASC

This is the code and data repository for the paper **Threshold-driven pruning with segmented maximum term weights for approximate cluster-based sparse retrieval** by Yifan Qiao, Parker Carlson, Shanxiu He, Yingrui Yang, and Tao Yang. 

The code is forked from and build upon [PISA](https://github.com/pisa-engine/pisa).

## Instructions

To build the project,

```
mkdir build
cd build
cmake ..
make -j reorder-docids lexicon compress_inverted_index create_wand_data evaluate_queries
```

To execute the code,

```
./evaluate_queries -e block_simdbp -i $COMPRESSED_INDEX_FILE -w $WAND_FILE -q $QUERY_FILE -k 1000 -a maxscore_boun
dsum --weighted -s quantized --documents $DOC_LEX_FILE --max-clusters $CLUSTER_SIZE --threads 0 --mu $MU_VAL —eta $ETA_VAL --subcluster_size_path $SUB_CLUSTER_SIZE_FILE
```

The required files can either be downloaded from [here](https://drive.google.com/drive/folders/1La2PL27fjYbjp8dtd0lV_G9HK4ZjMpF-?usp=sharing), or generated from the following commands:

```
./reorder-docids --collection $SPLADE_ORIGINAL_FILE --output $REORDERED_INDEX_PATH --from-mapping $CLUSTER_MAPPING_FILE

./lexicon build $DOC_MAP_FILE $DOC_LEX_FILE

./compress_inverted_index -c $REORDERED_INDEX_PATH -o $COMPRESSED_INDEX_FILE -e block_simdbp

./create_wand_data -c $REORDERED_INDEX_PATH -o $WAND_FILE -b 4000 -s quantized --document-clusters $CLUSTER_INFO_FILE
```

## Reference

If ASC helps you in your research, you can cite the following reference:
```
@inproceedings{qiao-etal-2024-threshold,
    title = "Threshold-driven Pruning with Segmented Maximum Term Weights for Approximate Cluster-based Sparse Retrieval",
    author = "Qiao, Yifan  and
      Carlson, Parker  and
      He, Shanxiu  and
      Yang, Yingrui  and
      Yang, Tao",
    editor = "Al-Onaizan, Yaser  and
      Bansal, Mohit  and
      Chen, Yun-Nung",
    booktitle = "Proceedings of the 2024 Conference on Empirical Methods in Natural Language Processing",
    month = nov,
    year = "2024",
    address = "Miami, Florida, USA",
    publisher = "Association for Computational Linguistics",
    url = "https://aclanthology.org/2024.emnlp-main.1101/",
    doi = "10.18653/v1/2024.emnlp-main.1101",
    pages = "19742--19757",
    abstract = "This paper revisits dynamic pruning through rank score thresholding in cluster-based sparse retrieval to skip the index partially at cluster and document levels during inference. It proposes a two-parameter pruning control scheme called ASC with a probabilistic guarantee on rank-safeness competitiveness. ASC uses cluster-level maximum weight segmentation to improve accuracy of rank score bound estimation and threshold-driven pruning, and is targeted for speeding up retrieval applications requiring high relevance competitiveness. The experiments with MS MARCO and BEIR show that ASC improves the accuracy and safeness of pruning for better relevance while delivering a low latency on a single-threaded CPU."
}
```
