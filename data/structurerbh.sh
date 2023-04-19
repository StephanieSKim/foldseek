#!/bin/sh -e
# reciprocal best hit workflow
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <sequenceDB> <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[   -f "$3.dbtype" ] && echo "$3.dbtype exists already!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";

A_DB="$1"
B_DB="$2"
RBH_RES="$3"
TMP_PATH="$4"

# search in both directions:
if [ ! -e "${TMP_PATH}/resAB.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" search "${A_DB}" "${B_DB}" "${TMP_PATH}/resAB" "${TMP_PATH}/tempAB" --min-seq-id 0 -c 0 --cov-mode 0 --max-rejected 2147483647 --max-accept 2147483647 -a 1 --add-self-matches 0 --tmscore-threshold 0 --tmalign-hit-order 0 --tmalign-fast 1 --db-load-mode 0 --threads 16 -v 3 --lddt-threshold 0 --sort-by-structure-bits 0 --sub-mat 'aa:3di.out,nucl:3di.out' --alignment-mode 3 --alignment-output-mode 0 --wrapped-scoring 0 -e 10 --min-aln-len 0 --seq-id-mode 0 --alt-ali 0 --max-seq-len 65535 --comp-bias-corr 0 --comp-bias-corr-scale 1 --pca substitution:1.100,context:1.400 --pcb substitution:4.100,context:5.800 --score-bias 0 --realign 0 --realign-score-bias -0.2 --realign-max-seqs 2147483647 --corr-score-weight 0 --gap-open aa:10,nucl:10 --gap-extend aa:1,nucl:1 --zdrop 40 --compressed 0 --seed-sub-mat 'aa:3di.out,nucl:3di.out' -s 9.5 -k 0 --k-score seq:2147483647,prof:2147483647 --alph-size aa:21,nucl:5 --max-seqs 1000 --split 0 --split-mode 2 --split-memory-limit 0 --diag-score 1 --exact-kmer-matching 0 --mask 0 --mask-prob 0.99995 --mask-lower-case 1 --min-ungapped-score 15 --spaced-kmer-mode 1 --exhaustive-search 0 --num-iterations 1 --alignment-type 2 --remove-tmp-files 0 --force-reuse 0 \
        || fail "search A vs. B died"
fi

# swap the direction of resBA::
if [ ! -e "${TMP_PATH}/resBA_swap.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" swapresults "${A_DB}" "${B_DB}" "${TMP_PATH}/tempAB/latest/pref" "${TMP_PATH}/prefBA" ${THREADS_COMP_PAR} -e inf \
        || fail "swap B prefilter A died"
fi

# structalign resBA:
if [ ! -e "${TMP_PATH}/resBA.dbtype" ]; then
    "$MMSEQS" structurealign "${B_DB}" "${A_DB}"  "${TMP_PATH}/prefBA" "${TMP_PATH}/resBA" --threads 16 -a \
       || fail "align died"
fi

#if [ ! -e "${TMP_PATH}/resBA.dbtype" ]; then
#    # shellcheck disable=SC2086
#    "$MMSEQS" search "${B_DB}" "${A_DB}" "${TMP_PATH}/resBA" "${TMP_PATH}/tempBA" ${SEARCH_B_A_PAR} \
#        || fail "search B vs. A died"
#fi


if [ ! -e "${TMP_PATH}/resAB_rescored.dbtype" ]; then
      # shellcheck disable=SC2086
    "$MMSEQS" rescorebacktrace "${A_DB}" "${B_DB}" "${TMP_PATH}/resAB" "${TMP_PATH}/resAB_rescored" ${THREADS_COMP_PAR} -a
fi
if [ ! -e "${TMP_PATH}/resBA_rescored.dbtype" ]; then
      # shellcheck disable=SC2086
  "$MMSEQS" rescorebacktrace "${B_DB}" "${A_DB}" "${TMP_PATH}/resBA" "${TMP_PATH}/resBA_rescored" ${THREADS_COMP_PAR} -a
fi

# extract a single best hit in A->B direction (used to take best bitscore for A):
if [ ! -e "${TMP_PATH}/resA_best_B.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" filterdb "${TMP_PATH}/resAB_rescored" "${TMP_PATH}/resA_best_B" --extract-lines 1 ${THREADS_COMP_PAR} \
        || fail "extract A best B died"
fi

# extract best hit(s) in B->A direction:
if [ ! -e "${TMP_PATH}/resB_best_A.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" filterdb "${TMP_PATH}/resBA_rescored" "${TMP_PATH}/resB_best_A" --beats-first --filter-column 2 --comparison-operator e ${THREADS_COMP_PAR} \
        || fail "extract B best A died"
fi

# swap the direction of resB_best_A:
if [ ! -e "${TMP_PATH}/resB_best_A_swap.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" swapresults "${B_DB}" "${A_DB}" "${TMP_PATH}/resB_best_A" "${TMP_PATH}/resB_best_A_swap" ${THREADS_COMP_PAR} -e 100000000 \
        || fail "swap B best A died"
fi

# merge the best results:
if [ ! -e "${TMP_PATH}/res_best_merged.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" mergedbs "${TMP_PATH}/resA_best_B" "${TMP_PATH}/res_best_merged" "${TMP_PATH}/resA_best_B" "${TMP_PATH}/resB_best_A_swap" ${VERB_COMP_PAR} \
        || fail "merge best hits died"
fi

# identify the RBH pairs and write them to a result db:
if [ ! -e "${RBH_RES}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" result2rbh "${TMP_PATH}/res_best_merged" "${RBH_RES}" ${THREADS_COMP_PAR} \
        || fail "result2rbh died"
fi

if [ -n "$REMOVE_TMP" ]; then
    rm -rf "${TMP_PATH}/tempAB"
    rm -rf "${TMP_PATH}/tempBA"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/resAB" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/resBA" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/resAB_rescored" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/resBA_rescored" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/resA_best_B" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/resB_best_A" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/res_best" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/resB_best_A_swap" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/res_best_merged_sorted" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/res_best_merged" ${VERBOSITY}
    rm -f "${TMP_PATH}/rbh.sh"
fi
