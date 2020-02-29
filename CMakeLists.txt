cmake_minimum_required(VERSION 3.2)

project(popscle VERSION 1.1.0)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;${CMAKE_FIND_LIBRARY_SUFFIXES}") # Prefer libz.a when both are available

find_package(ZLIB REQUIRED)
find_library(HTS_LIBRARY hts)
find_library(ZLIB_LIBRARY z)
find_package(Threads)

add_executable(popscle
        allele.cpp
        allele.h
        bam_ordered_reader.cpp
        bam_ordered_reader.h
        bcf_chunked_reader.cpp
        bcf_chunked_reader.h
        bcf_filter_arg.h
        bcf_filtered_reader.cpp
        bcf_filtered_reader.h
        bcf_variant_key.h
        bgzf.h
        cmd_cram_demuxlet.cpp
        cmd_cram_digital_pileup.cpp
        cmd_cram_freemuxlet.cpp
        cmd_cram_freemux2.cpp
        cmd_plp_make_dge_matrix.cpp
        commands.cpp
        commands.h
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
        )

target_include_directories(popscle PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)

target_link_libraries(popscle ${HTS_LIBRARY} ${ZLIB_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS popscle DESTINATION bin)
