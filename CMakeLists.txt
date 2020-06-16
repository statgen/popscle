cmake_minimum_required(VERSION 2.8.4)
project(popscle)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")


find_path(HTS_INCLUDE_DIRS htslib/bgzf.h HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT HTS_INCLUDE_DIRS )
    message(FATAL_ERROR "libhts HTS_INCLUDE_DIRS not found")
endif()

message(STATUS ${HTS_INCLUDE_DIRS})

find_library(HTS_LIBRARIES hts HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT HTS_LIBRARIES)
    message(FATAL_ERROR "libhts HTS_LIBRARIES not found")
endif()

message(STATUS ${HTS_LIBRARIES})

include_directories(${HTS_INCLUDE_DIRS})

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
        cmd_cram_dsc_pileup.cpp
        cmd_cram_dsc_dump.cpp
        cmd_cram_dsc_dump2pileup.cpp	
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
	bcftools.h
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

find_library(ZLIB z HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT ZLIB)
    message(FATAL_ERROR "libz library not found")
endif()

find_library(CURLLIB curl HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CURLLIB)
    message(FATAL_ERROR "libcurl library not found")
endif()

find_package (BZip2)
if (NOT BZIP2_FOUND)
    message(FATAL_ERROR "libbz2 library not found")
else()
    include_directories(${BZIP2_INCLUDE_DIRS})
    set(BZIP2 ${BZIP2_LIBRARIES})
endif()

find_library(LZMA lzma HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT LZMA)
    message(WARNING "liblzma library not found, if you specified --disable-lzma when compiling libhts, please ignore this warning")
    set(LZMA "")
else()
    set(LZMA ${LZMA_LIBRARIES})
endif()

find_library(CRYPTOLIB crypto HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CRYPTOLIB)
    message(FATAL_ERROR "libcrypto library not found")
endif()

target_link_libraries(popscle ${HTS_LIBRARIES} ${ZLIB} ${BZIP2} ${LZMA} ${CURLLIB} ${CRYPTOLIB})

install(TARGETS popscle RUNTIME DESTINATION bin)
