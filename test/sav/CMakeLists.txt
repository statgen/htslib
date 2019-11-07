cmake_minimum_required(VERSION 3.2)
project(test-sav VERSION 1.0.0)

include_directories(SYSTEM ${CMAKE_CURRENT_LIST_DIR}/../..)
add_definitions(-DHTS_VERSION_TEXT="1.9" -DPLUGIN_EXT=".so")

add_executable(test-sav test_sav.c
        ../../bcf_sr_sort.c ../../bgzf.c ../../errmod.c ../../faidx.c ../../header.c ../../hfile.c ../../hfile_gcs.c ../../hfile_libcurl.c ../../hfile_net.c ../../hfile_s3.c ../../hfile_s3_write.c ../../hts.c ../../hts_os.c ../../kfunc.c ../../knetfile.c ../../kstring.c ../../md5.c ../../multipart.c ../../plugin.c ../../probaln.c ../../realn.c ../../regidx.c ../../region.c ../../sam.c ../../synced_bcf_reader.c ../../tbx.c ../../textutils.c ../../thread_pool.c ../../vcf.c ../../vcf_sweep.c ../../vcfutils.c
        ../../cram/cram_codecs.c ../../cram/cram_decode.c ../../cram/cram_encode.c ../../cram/cram_external.c ../../cram/cram_index.c ../../cram/cram_io.c ../../cram/cram_samtools.c ../../cram/cram_stats.c ../../cram/mFILE.c ../../cram/open_trace_file.c ../../cram/pooled_alloc.c ../../cram/rANS_static.c ../../cram/string_alloc.c
        )
target_link_libraries(test-sav z curl)