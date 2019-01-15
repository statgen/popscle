cmake_minimum_required(VERSION 2.8.4)
project(popscle)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")


include(ExternalProject)
ExternalProject_Add(HTSLIB
        PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
        STAMP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/build/stamp
        GIT_REPOSITORY https://github.com/samtools/htslib.git
        GIT_TAG "master"
        UPDATE_COMMAND ""
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/htslib
        CONFIGURE_COMMAND autoheader && autoconf && ${CMAKE_CURRENT_SOURCE_DIR}/htslib/configure --prefix=${CMAKE_CURRENT_SOURCE_DIR}/htslib --disable-lzma
        BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/htslib
        BUILD_COMMAND make && make test
        )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/htslib/include/)

add_executable(popscle
        allele.cpp
        allele.h
        bam_endian.h
        bam_ordered_reader.cpp
        bam_ordered_reader.h
        bcf_chunked_reader.cpp
        bcf_chunked_reader.h
        bcf_filter_arg.h
        bcf_filtered_reader.cpp
        bcf_filtered_reader.h
        bcf_variant_key.h
        bgzf.h
        bh_tsne.cpp
        bh_tsne.h
        cmd_cram_demuxlet.cpp
        cmd_cram_digital_pileup.cpp
        cmd_cram_freemuxlet.cpp
        cmd_cram_gtf_util.cpp
        cmd_cram_mux_pileup.cpp
        cmd_cram_simuxlet.cpp
        cmd_plp_find_mux_clust.cpp
        cmd_plp_make_dge_matrix.cpp
        cmd_sc_kallisto_count.cpp
        cmd_sc_map_stamps.cpp
        cmd_sc_multinom_em.cpp
        cmd_sc_multinom_gibbs.cpp
        commands.cpp
        commands.h
        config.h
        Constant.h
        cramore.cpp
        cramore.h
        discrete_log_helper.h
        dropseq.cpp
        dropseq.h
        Error.cpp
        Error.h
        filter.cpp
        filter.h
        genome_interval.cpp
        genome_interval.h
        genomeChunk.cpp
        genomeChunk.h
        genomeLoci.h
        gtf.cpp
        gtf.h
        gtf_interval_tree.h
        hts_utils.cpp
        hts_utils.h
        interval.cpp
        interval.h
        interval_tree.cpp
        interval_tree.h
        IO.cpp
        IO.h
        knetfile.h
        ksort.h
        kstring.h
        log_tool.cpp
        log_tool.h
        louvain.h
        params.cpp
        params.h
        pException.h
        PhredHelper.cpp
        PhredHelper.h
        reference_sequence.cpp
        reference_sequence.h
        sam_filtered_reader.cpp
        sam_filtered_reader.h
        sam_ordered_writer.cpp
        sam_ordered_writer.h
        sc_drop_seq.cpp
        sc_drop_seq.h
        sptree.cpp
        sptree.h
        tcf.h
        tsv_reader.cpp
        tsv_reader.h
        utils.cpp
        utils.h
        variant.cpp
        variant.h
        variant_manip.cpp
        variant_manip.h
        vntr.cpp
        vntr.h
        vptree.h
        )
add_dependencies(popscle HTSLIB)

find_library(ZLIB z HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT ZLIB)
    message(FATAL_ERROR "libz library not found")
endif()

find_library(CURLLIB curl HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CURLLIB)
    message(FATAL_ERROR "libcurl library not found")
endif()

find_library(BZ2LIB bz2 HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT BZ2LIB)
    message(FATAL_ERROR "libbz2 library not found")
endif()

find_library(CRYPTOLIB crypto HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CRYPTOLIB)
    message(FATAL_ERROR "libcrypto library not found")
endif()

target_link_libraries(popscle ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a ${ZLIB} ${BZ2LIB} ${CURLLIB} ${CRYPTOLIB})
