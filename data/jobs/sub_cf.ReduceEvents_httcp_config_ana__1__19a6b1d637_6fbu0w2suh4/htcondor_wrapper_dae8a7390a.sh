#!/usr/bin/env bash

# Wrapper script that is to be configured as htcondor's main executable file

htcondor_wrapper() {
    # helper to select the correct python executable
    _law_python() {
        command -v python &> /dev/null && python "$@" || python3 "$@"
    }

    #
    # detect variables
    #

    local shell_is_zsh="$( [ -z "${ZSH_VERSION}" ] && echo "false" || echo "true" )"
    local this_file="$( ${shell_is_zsh} && echo "${(%):-%x}" || echo "${BASH_SOURCE[0]}" )"
    local this_file_base="$( basename "${this_file}" )"

    # get the job number
    export LAW_HTCONDOR_JOB_NUMBER="${LAW_HTCONDOR_JOB_PROCESS}"
    if [ -z "${LAW_HTCONDOR_JOB_NUMBER}" ]; then
        >&2 echo "could not determine htcondor job number"
        return "1"
    fi
    # htcondor process numbers start at 0, law job numbers at 1, so increment
    ((LAW_HTCONDOR_JOB_NUMBER++))
    echo "running ${this_file_base} for job number ${LAW_HTCONDOR_JOB_NUMBER}"


    #
    # job argument definitons, depending on LAW_HTCONDOR_JOB_NUMBER
    #

    # definition
    local htcondor_job_arguments_map
    declare -A htcondor_job_arguments_map
    htcondor_job_arguments_map=(
        ['1']="columnflow.tasks.reduction ReduceEvents LS1sb2ctZmlsZT1OT19TVFIgLS1jbGVhci1sb2dzPUZhbHNlIC0tYW5hbHlzaXM9aHR0Y3AuY29uZmlnLmFuYWx5c2lzX2h0dGNwLmFuYWx5c2lzX2h0dGNwIC0tdmVyc2lvbj1qZXRfdmV0b19maXggLS1jb25maWc9cnVuM18yMDIyX3ByZUVFX211dGF1IC0tc2hpZnQ9bm9taW5hbCAtLWRhdGFzZXQ9aF9nZ2ZfaHR0X21tX2ZpbHRlcmVkIC0tY2FsaWJyYXRvcnM9bWFpbiAtLXNlbGVjdG9yPW1haW4gLS1zZWxlY3Rvci1zdGVwcz0nJyAtLWNoZWNrLWZpbml0ZS1vdXRwdXQ9RmFsc2UgLS1jaGVjay1vdmVybGFwcGluZy1pbnB1dHM9RmFsc2UgLS1jZi5DYWxpYnJhdGVFdmVudHMtd29ya2Zsb3c9aHRjb25kb3IgLS1jZi5DYWxpYnJhdGVFdmVudHMtdmVyc2lvbj1qZXRfdmV0b19maXggLS1jZi5TZWxlY3RFdmVudHMtd29ya2Zsb3c9aHRjb25kb3IgLS1jZi5TZWxlY3RFdmVudHMtdmVyc2lvbj1qZXRfdmV0b19maXggLS1jZi5NZXJnZVJlZHVjZWRFdmVudHMtd29ya2Zsb3c9aHRjb25kb3IgLS1jZi5NZXJnZVJlZHVjZWRFdmVudHMtdmVyc2lvbj1qZXRfdmV0b19maXggLS1jZi5NZXJnZVNlbGVjdGlvblN0YXRzLXZlcnNpb249amV0X3ZldG9fZml4IC0tY2YuUHJvdmlkZVJlZHVjZWRFdmVudHMtdmVyc2lvbj1qZXRfdmV0b19maXggLS1sb2NhbC1zY2hlZHVsZXI9VHJ1ZQ== MTA= 1 no LQ=="
        ['2']="columnflow.tasks.reduction ReduceEvents LS1sb2ctZmlsZT1OT19TVFIgLS1jbGVhci1sb2dzPUZhbHNlIC0tYW5hbHlzaXM9aHR0Y3AuY29uZmlnLmFuYWx5c2lzX2h0dGNwLmFuYWx5c2lzX2h0dGNwIC0tdmVyc2lvbj1qZXRfdmV0b19maXggLS1jb25maWc9cnVuM18yMDIyX3ByZUVFX211dGF1IC0tc2hpZnQ9bm9taW5hbCAtLWRhdGFzZXQ9aF9nZ2ZfaHR0X21tX2ZpbHRlcmVkIC0tY2FsaWJyYXRvcnM9bWFpbiAtLXNlbGVjdG9yPW1haW4gLS1zZWxlY3Rvci1zdGVwcz0nJyAtLWNoZWNrLWZpbml0ZS1vdXRwdXQ9RmFsc2UgLS1jaGVjay1vdmVybGFwcGluZy1pbnB1dHM9RmFsc2UgLS1jZi5DYWxpYnJhdGVFdmVudHMtd29ya2Zsb3c9aHRjb25kb3IgLS1jZi5DYWxpYnJhdGVFdmVudHMtdmVyc2lvbj1qZXRfdmV0b19maXggLS1jZi5TZWxlY3RFdmVudHMtd29ya2Zsb3c9aHRjb25kb3IgLS1jZi5TZWxlY3RFdmVudHMtdmVyc2lvbj1qZXRfdmV0b19maXggLS1jZi5NZXJnZVJlZHVjZWRFdmVudHMtd29ya2Zsb3c9aHRjb25kb3IgLS1jZi5NZXJnZVJlZHVjZWRFdmVudHMtdmVyc2lvbj1qZXRfdmV0b19maXggLS1jZi5NZXJnZVNlbGVjdGlvblN0YXRzLXZlcnNpb249amV0X3ZldG9fZml4IC0tY2YuUHJvdmlkZVJlZHVjZWRFdmVudHMtdmVyc2lvbj1qZXRfdmV0b19maXggLS1sb2NhbC1zY2hlZHVsZXI9VHJ1ZQ== MTE= 1 no LQ=="
        ['3']="columnflow.tasks.reduction ReduceEvents LS1sb2ctZmlsZT1OT19TVFIgLS1jbGVhci1sb2dzPUZhbHNlIC0tYW5hbHlzaXM9aHR0Y3AuY29uZmlnLmFuYWx5c2lzX2h0dGNwLmFuYWx5c2lzX2h0dGNwIC0tdmVyc2lvbj1qZXRfdmV0b19maXggLS1jb25maWc9cnVuM18yMDIyX3ByZUVFX211dGF1IC0tc2hpZnQ9bm9taW5hbCAtLWRhdGFzZXQ9aF9nZ2ZfaHR0X21tX2ZpbHRlcmVkIC0tY2FsaWJyYXRvcnM9bWFpbiAtLXNlbGVjdG9yPW1haW4gLS1zZWxlY3Rvci1zdGVwcz0nJyAtLWNoZWNrLWZpbml0ZS1vdXRwdXQ9RmFsc2UgLS1jaGVjay1vdmVybGFwcGluZy1pbnB1dHM9RmFsc2UgLS1jZi5DYWxpYnJhdGVFdmVudHMtd29ya2Zsb3c9aHRjb25kb3IgLS1jZi5DYWxpYnJhdGVFdmVudHMtdmVyc2lvbj1qZXRfdmV0b19maXggLS1jZi5TZWxlY3RFdmVudHMtd29ya2Zsb3c9aHRjb25kb3IgLS1jZi5TZWxlY3RFdmVudHMtdmVyc2lvbj1qZXRfdmV0b19maXggLS1jZi5NZXJnZVJlZHVjZWRFdmVudHMtd29ya2Zsb3c9aHRjb25kb3IgLS1jZi5NZXJnZVJlZHVjZWRFdmVudHMtdmVyc2lvbj1qZXRfdmV0b19maXggLS1jZi5NZXJnZVNlbGVjdGlvblN0YXRzLXZlcnNpb249amV0X3ZldG9fZml4IC0tY2YuUHJvdmlkZVJlZHVjZWRFdmVudHMtdmVyc2lvbj1qZXRfdmV0b19maXggLS1sb2NhbC1zY2hlZHVsZXI9VHJ1ZQ== MTI= 1 no LQ=="
    )

    # pick
    local htcondor_job_arguments="${htcondor_job_arguments_map[${LAW_HTCONDOR_JOB_NUMBER}]}"
    if [ -z "${htcondor_job_arguments}" ]; then
        >&2 echo "empty htcondor job arguments for LAW_HTCONDOR_JOB_NUMBER ${LAW_HTCONDOR_JOB_NUMBER}"
        return "3"
    fi


    #
    # variable rendering
    #

    # check variables
    local render_variables="eyJsYXdfY29uZmlnX2ZpbGUiOiAiJENGX1JFUE9fQkFTRS9sYXcuY2ZnIiwgImNmX3JlcG9fdXJpcyI6ICJyb290Oi8vZW9zdXNlci5jZXJuLmNoLy9lb3MvcHJvamVjdC9pL2lwaGN0YXUvcHVibGljL213aXR0L0NQaW5IVG9UYXVUYXVPdXRwdXQvaGFtYm91cmcvY2Zfc3RvcmUvYW5hbHlzaXNfaHR0Y3AvY2YuQnVuZGxlUmVwbyIsICJjZl9yZXBvX3BhdHRlcm4iOiAiQ1BpbkhUb1RhdVRhdS4zM2E1OTU1Y2YwYjNjZmYwNjhhZWUyOTk5Y2I3MzIxYjA3M2NhODZjLlteXFwuXSsudGd6IiwgImNmX3NvZnR3YXJlX3VyaXMiOiAicm9vdDovL2Vvc3VzZXIuY2Vybi5jaC8vZW9zL3Byb2plY3QvaS9pcGhjdGF1L3B1YmxpYy9td2l0dC9DUGluSFRvVGF1VGF1T3V0cHV0L2hhbWJvdXJnL2NmX3N0b3JlL2FuYWx5c2lzX2h0dGNwL2NmLkJ1bmRsZVNvZnR3YXJlIiwgImNmX3NvZnR3YXJlX3BhdHRlcm4iOiAic29mdHdhcmUuW15cXC5dKy50Z3oiLCAiY2ZfYmFzaF9zYW5kYm94X3VyaXMiOiAiXCJyb290Oi8vZW9zdXNlci5jZXJuLmNoLy9lb3MvcHJvamVjdC9pL2lwaGN0YXUvcHVibGljL213aXR0L0NQaW5IVG9UYXVUYXVPdXRwdXQvaGFtYm91cmcvY2Zfc3RvcmUvYW5hbHlzaXNfaHR0Y3AvY2YuQnVuZGxlQmFzaFNhbmRib3hcIiBcInJvb3Q6Ly9lb3N1c2VyLmNlcm4uY2gvL2Vvcy9wcm9qZWN0L2kvaXBoY3RhdS9wdWJsaWMvbXdpdHQvQ1BpbkhUb1RhdVRhdU91dHB1dC9oYW1ib3VyZy9jZl9zdG9yZS9hbmFseXNpc19odHRjcC9jZi5CdW5kbGVCYXNoU2FuZGJveFwiIiwgImNmX2Jhc2hfc2FuZGJveF9wYXR0ZXJucyI6ICJcImNmXzVkZTcyMTNjLjFkZmE0YWYyYTkuW15cXC5dKy50Z3pcIiBcInZlbnZfY29sdW1uYXJfNGVmNDU2NzQuNzI5MjUzNGRmZi5bXlxcLl0rLnRnelwiIiwgImNmX2Jhc2hfc2FuZGJveF9uYW1lcyI6ICJcImNmXCIgXCJ2ZW52X2NvbHVtbmFyXCIiLCAiY2ZfY21zc3dfc2FuZGJveF91cmlzIjogIlwicm9vdDovL2Vvc3VzZXIuY2Vybi5jaC8vZW9zL3Byb2plY3QvaS9pcGhjdGF1L3B1YmxpYy9td2l0dC9DUGluSFRvVGF1VGF1T3V0cHV0L2hhbWJvdXJnL2NmX3N0b3JlL2FuYWx5c2lzX2h0dGNwL2NmLkJ1bmRsZUNNU1NXU2FuZGJveFwiIiwgImNmX2Ntc3N3X3NhbmRib3hfcGF0dGVybnMiOiAiXCJjbXNzd19kZWZhdWx0Xzc1MGY2MWQ1X0NNU1NXXzE0XzFfMF9wcmU0LmRhMzlhM2VlNWU2YjRiMGQzMjU1YmZlZjk1NjAxODkwYWZkODA3MDkuW15cXC5dKy50Z3pcIiIsICJjZl9jbXNzd19zYW5kYm94X25hbWVzIjogIlwiY21zc3dfZGVmYXVsdFwiIiwgImNmX2Jvb3RzdHJhcF9uYW1lIjogImh0Y29uZG9yX3N0YW5kYWxvbmUiLCAiY2ZfaHRjb25kb3JfZmxhdm9yIjogImNlcm4iLCAiY2ZfcHJlX3NldHVwX2NvbW1hbmQiOiAiIiwgImNmX3Bvc3Rfc2V0dXBfY29tbWFuZCI6ICIiLCAiY2ZfcmVtb3RlX2xjZ19zZXR1cCI6ICIvY3ZtZnMvZ3JpZC5jZXJuLmNoL2FsbWE5LXVpLXRlc3QvZXRjL3Byb2ZpbGUuZC9zZXR1cC1hbG1hOS10ZXN0LnNoIiwgImNmX3JlbW90ZV9sY2dfc2V0dXBfZm9yY2UiOiAiIiwgImNmX2Jhc2UiOiAiL2Fmcy9jZXJuLmNoL3VzZXIvbS9td2l0dC9wdWJsaWMvQ1BpbkhUb1RhdVRhdS9tb2R1bGVzL2NvbHVtbmZsb3ciLCAiY2ZfcmVwb19iYXNlIjogIi9hZnMvY2Vybi5jaC91c2VyL20vbXdpdHQvcHVibGljL0NQaW5IVG9UYXVUYXUiLCAiY2ZfY2Vybl91c2VyIjogIm13aXR0IiwgImNmX3N0b3JlX25hbWUiOiAiY2Zfc3RvcmUiLCAiY2Zfc3RvcmVfbG9jYWwiOiAiL2Fmcy9jZXJuLmNoL3VzZXIvbS9td2l0dC9wdWJsaWMvQ1BpbkhUb1RhdVRhdS9kYXRhL2NmX3N0b3JlIiwgImNmX2xvY2FsX3NjaGVkdWxlciI6ICJ0cnVlIiwgImpvYl9maWxlIjogImxhd19qb2JfZDA4NTk5YmMwOC5zaCIsICJleGVjdXRhYmxlX2ZpbGUiOiAiaHRjb25kb3Jfd3JhcHBlcl9kYWU4YTczOTBhLnNoIiwgImJvb3RzdHJhcF9maWxlIjogInJlbW90ZV9ib290c3RyYXBfZDA3ZTUwYmE5NS5zaCIsICJ2b21zcHJveHlfZmlsZSI6ICJ4NTA5dXBfdTE1MTc5MF9hOGYxMmZkMWY1IiwgIndsY2dfdG9vbHMiOiAibGF3X3dsY2dfdG9vbHNfNDAxZDliOTA1Ny5zaCIsICJpbnB1dF9maWxlcyI6ICJsYXdfam9iX2QwODU5OWJjMDguc2ggaHRjb25kb3Jfd3JhcHBlcl9kYWU4YTczOTBhLnNoIHJlbW90ZV9ib290c3RyYXBfZDA3ZTUwYmE5NS5zaCB4NTA5dXBfdTE1MTc5MF9hOGYxMmZkMWY1IGxhd193bGNnX3Rvb2xzXzQwMWQ5YjkwNTcuc2giLCAiaW5wdXRfZmlsZXNfcmVuZGVyIjogImxhd19qb2JfZDA4NTk5YmMwOC5zaCByZW1vdGVfYm9vdHN0cmFwX2QwN2U1MGJhOTUuc2giLCAibG9nX2ZpbGUiOiAic3RkYWxsJChsYXdfam9iX3Bvc3RmaXgpLnR4dCIsICJodGNvbmRvcl9qb2JfYXJndW1lbnRzX21hcCI6ICJbJzEnXT1cImNvbHVtbmZsb3cudGFza3MucmVkdWN0aW9uIFJlZHVjZUV2ZW50cyBMUzFzYjJjdFptbHNaVDFPVDE5VFZGSWdMUzFqYkdWaGNpMXNiMmR6UFVaaGJITmxJQzB0WVc1aGJIbHphWE05YUhSMFkzQXVZMjl1Wm1sbkxtRnVZV3g1YzJselgyaDBkR053TG1GdVlXeDVjMmx6WDJoMGRHTndJQzB0ZG1WeWMybHZiajFxWlhSZmRtVjBiMTltYVhnZ0xTMWpiMjVtYVdjOWNuVnVNMTh5TURJeVgzQnlaVVZGWDIxMWRHRjFJQzB0YzJocFpuUTlibTl0YVc1aGJDQXRMV1JoZEdGelpYUTlhRjluWjJaZmFIUjBYMjF0WDJacGJIUmxjbVZrSUMwdFkyRnNhV0p5WVhSdmNuTTliV0ZwYmlBdExYTmxiR1ZqZEc5eVBXMWhhVzRnTFMxelpXeGxZM1J2Y2kxemRHVndjejBuSnlBdExXTm9aV05yTFdacGJtbDBaUzF2ZFhSd2RYUTlSbUZzYzJVZ0xTMWphR1ZqYXkxdmRtVnliR0Z3Y0dsdVp5MXBibkIxZEhNOVJtRnNjMlVnTFMxalppNURZV3hwWW5KaGRHVkZkbVZ1ZEhNdGQyOXlhMlpzYjNjOWFIUmpiMjVrYjNJZ0xTMWpaaTVEWVd4cFluSmhkR1ZGZG1WdWRITXRkbVZ5YzJsdmJqMXFaWFJmZG1WMGIxOW1hWGdnTFMxalppNVRaV3hsWTNSRmRtVnVkSE10ZDI5eWEyWnNiM2M5YUhSamIyNWtiM0lnTFMxalppNVRaV3hsWTNSRmRtVnVkSE10ZG1WeWMybHZiajFxWlhSZmRtVjBiMTltYVhnZ0xTMWpaaTVOWlhKblpWSmxaSFZqWldSRmRtVnVkSE10ZDI5eWEyWnNiM2M5YUhSamIyNWtiM0lnTFMxalppNU5aWEpuWlZKbFpIVmpaV1JGZG1WdWRITXRkbVZ5YzJsdmJqMXFaWFJmZG1WMGIxOW1hWGdnTFMxalppNU5aWEpuWlZObGJHVmpkR2x2YmxOMFlYUnpMWFpsY25OcGIyNDlhbVYwWDNabGRHOWZabWw0SUMwdFkyWXVVSEp2ZG1sa1pWSmxaSFZqWldSRmRtVnVkSE10ZG1WeWMybHZiajFxWlhSZmRtVjBiMTltYVhnZ0xTMXNiMk5oYkMxelkyaGxaSFZzWlhJOVZISjFaUT09IE1UQT0gMSBubyBMUT09XCJcbiAgICAgICAgWycyJ109XCJjb2x1bW5mbG93LnRhc2tzLnJlZHVjdGlvbiBSZWR1Y2VFdmVudHMgTFMxc2IyY3RabWxzWlQxT1QxOVRWRklnTFMxamJHVmhjaTFzYjJkelBVWmhiSE5sSUMwdFlXNWhiSGx6YVhNOWFIUjBZM0F1WTI5dVptbG5MbUZ1WVd4NWMybHpYMmgwZEdOd0xtRnVZV3g1YzJselgyaDBkR053SUMwdGRtVnljMmx2YmoxcVpYUmZkbVYwYjE5bWFYZ2dMUzFqYjI1bWFXYzljblZ1TTE4eU1ESXlYM0J5WlVWRlgyMTFkR0YxSUMwdGMyaHBablE5Ym05dGFXNWhiQ0F0TFdSaGRHRnpaWFE5YUY5bloyWmZhSFIwWDIxdFgyWnBiSFJsY21Wa0lDMHRZMkZzYVdKeVlYUnZjbk05YldGcGJpQXRMWE5sYkdWamRHOXlQVzFoYVc0Z0xTMXpaV3hsWTNSdmNpMXpkR1Z3Y3owbkp5QXRMV05vWldOckxXWnBibWwwWlMxdmRYUndkWFE5Um1Gc2MyVWdMUzFqYUdWamF5MXZkbVZ5YkdGd2NHbHVaeTFwYm5CMWRITTlSbUZzYzJVZ0xTMWpaaTVEWVd4cFluSmhkR1ZGZG1WdWRITXRkMjl5YTJac2IzYzlhSFJqYjI1a2IzSWdMUzFqWmk1RFlXeHBZbkpoZEdWRmRtVnVkSE10ZG1WeWMybHZiajFxWlhSZmRtVjBiMTltYVhnZ0xTMWpaaTVUWld4bFkzUkZkbVZ1ZEhNdGQyOXlhMlpzYjNjOWFIUmpiMjVrYjNJZ0xTMWpaaTVUWld4bFkzUkZkbVZ1ZEhNdGRtVnljMmx2YmoxcVpYUmZkbVYwYjE5bWFYZ2dMUzFqWmk1TlpYSm5aVkpsWkhWalpXUkZkbVZ1ZEhNdGQyOXlhMlpzYjNjOWFIUmpiMjVrYjNJZ0xTMWpaaTVOWlhKblpWSmxaSFZqWldSRmRtVnVkSE10ZG1WeWMybHZiajFxWlhSZmRtVjBiMTltYVhnZ0xTMWpaaTVOWlhKblpWTmxiR1ZqZEdsdmJsTjBZWFJ6TFhabGNuTnBiMjQ5YW1WMFgzWmxkRzlmWm1sNElDMHRZMll1VUhKdmRtbGtaVkpsWkhWalpXUkZkbVZ1ZEhNdGRtVnljMmx2YmoxcVpYUmZkbVYwYjE5bWFYZ2dMUzFzYjJOaGJDMXpZMmhsWkhWc1pYSTlWSEoxWlE9PSBNVEU9IDEgbm8gTFE9PVwiXG4gICAgICAgIFsnMyddPVwiY29sdW1uZmxvdy50YXNrcy5yZWR1Y3Rpb24gUmVkdWNlRXZlbnRzIExTMXNiMmN0Wm1sc1pUMU9UMTlUVkZJZ0xTMWpiR1ZoY2kxc2IyZHpQVVpoYkhObElDMHRZVzVoYkhsemFYTTlhSFIwWTNBdVkyOXVabWxuTG1GdVlXeDVjMmx6WDJoMGRHTndMbUZ1WVd4NWMybHpYMmgwZEdOd0lDMHRkbVZ5YzJsdmJqMXFaWFJmZG1WMGIxOW1hWGdnTFMxamIyNW1hV2M5Y25WdU0xOHlNREl5WDNCeVpVVkZYMjExZEdGMUlDMHRjMmhwWm5ROWJtOXRhVzVoYkNBdExXUmhkR0Z6WlhROWFGOW5aMlpmYUhSMFgyMXRYMlpwYkhSbGNtVmtJQzB0WTJGc2FXSnlZWFJ2Y25NOWJXRnBiaUF0TFhObGJHVmpkRzl5UFcxaGFXNGdMUzF6Wld4bFkzUnZjaTF6ZEdWd2N6MG5KeUF0TFdOb1pXTnJMV1pwYm1sMFpTMXZkWFJ3ZFhROVJtRnNjMlVnTFMxamFHVmpheTF2ZG1WeWJHRndjR2x1WnkxcGJuQjFkSE05Um1Gc2MyVWdMUzFqWmk1RFlXeHBZbkpoZEdWRmRtVnVkSE10ZDI5eWEyWnNiM2M5YUhSamIyNWtiM0lnTFMxalppNURZV3hwWW5KaGRHVkZkbVZ1ZEhNdGRtVnljMmx2YmoxcVpYUmZkbVYwYjE5bWFYZ2dMUzFqWmk1VFpXeGxZM1JGZG1WdWRITXRkMjl5YTJac2IzYzlhSFJqYjI1a2IzSWdMUzFqWmk1VFpXeGxZM1JGZG1WdWRITXRkbVZ5YzJsdmJqMXFaWFJmZG1WMGIxOW1hWGdnTFMxalppNU5aWEpuWlZKbFpIVmpaV1JGZG1WdWRITXRkMjl5YTJac2IzYzlhSFJqYjI1a2IzSWdMUzFqWmk1TlpYSm5aVkpsWkhWalpXUkZkbVZ1ZEhNdGRtVnljMmx2YmoxcVpYUmZkbVYwYjE5bWFYZ2dMUzFqWmk1TlpYSm5aVk5sYkdWamRHbHZibE4wWVhSekxYWmxjbk5wYjI0OWFtVjBYM1psZEc5ZlptbDRJQzB0WTJZdVVISnZkbWxrWlZKbFpIVmpaV1JGZG1WdWRITXRkbVZ5YzJsdmJqMXFaWFJmZG1WMGIxOW1hWGdnTFMxc2IyTmhiQzF6WTJobFpIVnNaWEk5VkhKMVpRPT0gTVRJPSAxIG5vIExRPT1cIiJ9"
    if [ -z "${render_variables}" ]; then
        >&2 echo "empty render variables"
        return "4"
    fi

    # decode
    render_variables="$( echo "${render_variables}" | base64 --decode )"

    # check files to render
    local input_files_render=( law_job_d08599bc08.sh remote_bootstrap_d07e50ba95.sh )
    if [ "${#input_files_render[@]}" == "0" ]; then
        >&2 echo "received empty input files for rendering for LAW_HTCONDOR_JOB_NUMBER ${LAW_HTCONDOR_JOB_NUMBER}"
        return "5"
    fi

    # render files
    local input_file_render
    for input_file_render in ${input_files_render[@]}; do
        # skip if the file refers to _this_ one
        local input_file_render_base="$( basename "${input_file_render}" )"
        [ "${input_file_render_base}" = "${this_file_base}" ] && continue
        # render
        echo "render ${input_file_render}"
        cat > _render.py << EOT
import re
repl = ${render_variables}
repl['input_files_render'] = ''
repl['file_postfix'] = '${file_postfix}' or repl.get('file_postfix', '')
repl['log_file'] = ''
content = open('${input_file_render}', 'r').read()
content = re.sub(r'\{\{(\w+)\}\}', lambda m: repl.get(m.group(1), ''), content)
open('${input_file_render_base}', 'w').write(content)
EOT
        _law_python _render.py
        local render_ret="$?"
        rm -f _render.py
        # handle rendering errors
        if [ "${render_ret}" != "0" ]; then
            >&2 echo "input file rendering failed with code ${render_ret}"
            return "6"
        fi
    done


    #
    # run the actual job file
    #

    # check the job file
    local job_file="law_job_d08599bc08.sh"
    if [ ! -f "${job_file}" ]; then
        >&2 echo "job file '${job_file}' does not exist"
        return "7"
    fi

    # helper to print a banner
    banner() {
        local msg="$1"

        echo
        echo "================================================================================"
        echo "=== ${msg}"
        echo "================================================================================"
        echo
    }

    # debugging: print its contents
    # echo "=== content of job file '${job_file}'"
    # echo
    # cat "${job_file}"
    # echo
    # echo "=== end of job file content"

    # run it
    banner "Start of law job"

    local job_ret
    bash "${job_file}" ${htcondor_job_arguments}
    job_ret="$?"

    banner "End of law job"

    return "${job_ret}"
}

action() {
    # arguments: file_postfix, log_file
    local file_postfix="$1"
    local log_file="$2"

    # create log directory
    if [ ! -z "${log_file}" ]; then
        local log_dir="$( dirname "${log_file}" )"
        [ ! -d "${log_dir}" ] && mkdir -p "${log_dir}"
    fi

    # run the wrapper function
    if [ -z "${log_file}" ]; then
        htcondor_wrapper "$@"
    elif command -v tee &> /dev/null; then
        set -o pipefail
        echo "---" >> "${log_file}"
        htcondor_wrapper "$@" 2>&1 | tee -a "${log_file}"
    else
        echo "---" >> "${log_file}"
        htcondor_wrapper "$@" &>> "${log_file}"
    fi
}

action "$@"
