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
        ['1']="columnflow.tasks.selection SelectEvents LS1sb2ctZmlsZT1OT19TVFIgLS1jbGVhci1sb2dzPUZhbHNlIC0tYW5hbHlzaXM9aHR0Y3AuY29uZmlnLmFuYWx5c2lzX2h0dGNwLmFuYWx5c2lzX2h0dGNwIC0tdmVyc2lvbj1mcmFtZXdvcmtfc3luYyAtLWNvbmZpZz1ydW4zXzIwMjJfcHJlRUVfbXV0YXUgLS1zaGlmdD1ub21pbmFsIC0tZGF0YXNldD1keV9sZXBfbWFkZ3JhcGggLS1jYWxpYnJhdG9ycz1tYWluIC0tc2VsZWN0b3I9bWFpbiAtLWNoZWNrLWZpbml0ZS1vdXRwdXQ9RmFsc2UgLS1jaGVjay1vdmVybGFwcGluZy1pbnB1dHM9RmFsc2UgLS1jZi5DYWxpYnJhdGVFdmVudHMtd29ya2Zsb3c9aHRjb25kb3IgLS1jZi5DYWxpYnJhdGVFdmVudHMtdmVyc2lvbj1mcmFtZXdvcmtfc3luYyAtLWNmLlJlZHVjZUV2ZW50cy13b3JrZmxvdz1odGNvbmRvciAtLWNmLlJlZHVjZUV2ZW50cy12ZXJzaW9uPWZyYW1ld29ya19zeW5jIC0tY2YuTWVyZ2VSZWR1Y2VkRXZlbnRzLXdvcmtmbG93PWh0Y29uZG9yIC0tY2YuTWVyZ2VSZWR1Y2VkRXZlbnRzLXZlcnNpb249ZnJhbWV3b3JrX3N5bmMgLS1jZi5NZXJnZVNlbGVjdGlvblN0YXRzLXZlcnNpb249ZnJhbWV3b3JrX3N5bmMgLS1jZi5Qcm92aWRlUmVkdWNlZEV2ZW50cy12ZXJzaW9uPWZyYW1ld29ya19zeW5jIC0tbG9jYWwtc2NoZWR1bGVyPVRydWU= Nw== 1 no LQ=="
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
    local render_variables="eyJsYXdfY29uZmlnX2ZpbGUiOiAiJENGX1JFUE9fQkFTRS9sYXcuY2ZnIiwgImNmX3JlcG9fdXJpcyI6ICJyb290Oi8vZW9zdXNlci5jZXJuLmNoLy9lb3MvcHJvamVjdC9pL2lwaGN0YXUvcHVibGljL213aXR0L0NQaW5IVG9UYXVUYXVPdXRwdXQvaGFtYm91cmcvY2Zfc3RvcmUvYW5hbHlzaXNfaHR0Y3AvY2YuQnVuZGxlUmVwbyIsICJjZl9yZXBvX3BhdHRlcm4iOiAiQ1BpbkhUb1RhdVRhdS5iOTY2MzUyNjc1YjVkNDk5OGI2OTFiMDdmNGY4MTQzZjliZmFiYmQ5LlteXFwuXSsudGd6IiwgImNmX3NvZnR3YXJlX3VyaXMiOiAicm9vdDovL2Vvc3VzZXIuY2Vybi5jaC8vZW9zL3Byb2plY3QvaS9pcGhjdGF1L3B1YmxpYy9td2l0dC9DUGluSFRvVGF1VGF1T3V0cHV0L2hhbWJvdXJnL2NmX3N0b3JlL2FuYWx5c2lzX2h0dGNwL2NmLkJ1bmRsZVNvZnR3YXJlIiwgImNmX3NvZnR3YXJlX3BhdHRlcm4iOiAic29mdHdhcmUuW15cXC5dKy50Z3oiLCAiY2ZfYmFzaF9zYW5kYm94X3VyaXMiOiAiXCJyb290Oi8vZW9zdXNlci5jZXJuLmNoLy9lb3MvcHJvamVjdC9pL2lwaGN0YXUvcHVibGljL213aXR0L0NQaW5IVG9UYXVUYXVPdXRwdXQvaGFtYm91cmcvY2Zfc3RvcmUvYW5hbHlzaXNfaHR0Y3AvY2YuQnVuZGxlQmFzaFNhbmRib3hcIiBcInJvb3Q6Ly9lb3N1c2VyLmNlcm4uY2gvL2Vvcy9wcm9qZWN0L2kvaXBoY3RhdS9wdWJsaWMvbXdpdHQvQ1BpbkhUb1RhdVRhdU91dHB1dC9oYW1ib3VyZy9jZl9zdG9yZS9hbmFseXNpc19odHRjcC9jZi5CdW5kbGVCYXNoU2FuZGJveFwiIFwicm9vdDovL2Vvc3VzZXIuY2Vybi5jaC8vZW9zL3Byb2plY3QvaS9pcGhjdGF1L3B1YmxpYy9td2l0dC9DUGluSFRvVGF1VGF1T3V0cHV0L2hhbWJvdXJnL2NmX3N0b3JlL2FuYWx5c2lzX2h0dGNwL2NmLkJ1bmRsZUJhc2hTYW5kYm94XCIiLCAiY2ZfYmFzaF9zYW5kYm94X3BhdHRlcm5zIjogIlwiY2ZfNWRlNzIxM2MuOTAwMzY1NzI0ZS5bXlxcLl0rLnRnelwiIFwidmVudl9jb2x1bW5hcl80ZWY0NTY3NC44MjZiNjIwNTk4LlteXFwuXSsudGd6XCIgXCJ2ZW52X2NvbHVtbmFyX3hnYl9mMDUxZjcyZC5jNzA3NDViNzA5LlteXFwuXSsudGd6XCIiLCAiY2ZfYmFzaF9zYW5kYm94X25hbWVzIjogIlwiY2ZcIiBcInZlbnZfY29sdW1uYXJcIiBcInZlbnZfY29sdW1uYXJfeGdiXCIiLCAiY2ZfYm9vdHN0cmFwX25hbWUiOiAiaHRjb25kb3Jfc3RhbmRhbG9uZSIsICJjZl9odGNvbmRvcl9mbGF2b3IiOiAiY2VybiIsICJjZl9wcmVfc2V0dXBfY29tbWFuZCI6ICIiLCAiY2ZfcG9zdF9zZXR1cF9jb21tYW5kIjogIiIsICJjZl9yZW1vdGVfbGNnX3NldHVwIjogIi9jdm1mcy9ncmlkLmNlcm4uY2gvYWxtYTktdWktdGVzdC9ldGMvcHJvZmlsZS5kL3NldHVwLWFsbWE5LXRlc3Quc2giLCAiY2ZfcmVtb3RlX2xjZ19zZXR1cF9mb3JjZSI6ICIiLCAiY2ZfYmFzZSI6ICIvYWZzL2Nlcm4uY2gvdXNlci9tL213aXR0L21hc3RyaGgvQ1BpbkhUb1RhdVRhdS9tb2R1bGVzL2NvbHVtbmZsb3ciLCAiY2ZfcmVwb19iYXNlIjogIi9hZnMvY2Vybi5jaC91c2VyL20vbXdpdHQvbWFzdHJoaC9DUGluSFRvVGF1VGF1IiwgImNmX2Nlcm5fdXNlciI6ICJtd2l0dCIsICJjZl9zdG9yZV9uYW1lIjogImNmX3N0b3JlIiwgImNmX3N0b3JlX2xvY2FsIjogIi9hZnMvY2Vybi5jaC91c2VyL20vbXdpdHQvbWFzdHJoaC9DUGluSFRvVGF1VGF1L2RhdGEvY2Zfc3RvcmUiLCAiY2ZfbG9jYWxfc2NoZWR1bGVyIjogInRydWUiLCAiam9iX2ZpbGUiOiAibGF3X2pvYl82MGFiZGViNDg3LnNoIiwgImV4ZWN1dGFibGVfZmlsZSI6ICJodGNvbmRvcl93cmFwcGVyX2I2Y2NkMjIwNzlfOS5zaCIsICJib290c3RyYXBfZmlsZSI6ICJyZW1vdGVfYm9vdHN0cmFwXzExZDFiNjVhOTguc2giLCAidm9tc3Byb3h5X2ZpbGUiOiAieDUwOXVwX3UxNTE3OTBfYThmMTJmZDFmNSIsICJ3bGNnX3Rvb2xzIjogImxhd193bGNnX3Rvb2xzX2MxNGNhZTBhZmYuc2giLCAiaW5wdXRfZmlsZXMiOiAibGF3X2pvYl82MGFiZGViNDg3LnNoIGh0Y29uZG9yX3dyYXBwZXJfYjZjY2QyMjA3OV85LnNoIHJlbW90ZV9ib290c3RyYXBfMTFkMWI2NWE5OC5zaCB4NTA5dXBfdTE1MTc5MF9hOGYxMmZkMWY1IGxhd193bGNnX3Rvb2xzX2MxNGNhZTBhZmYuc2giLCAiaW5wdXRfZmlsZXNfcmVuZGVyIjogImxhd19qb2JfNjBhYmRlYjQ4Ny5zaCByZW1vdGVfYm9vdHN0cmFwXzExZDFiNjVhOTguc2giLCAibG9nX2ZpbGUiOiAic3RkYWxsJChsYXdfam9iX3Bvc3RmaXgpLnR4dCIsICJodGNvbmRvcl9qb2JfYXJndW1lbnRzX21hcCI6ICJbJzEnXT1cImNvbHVtbmZsb3cudGFza3Muc2VsZWN0aW9uIFNlbGVjdEV2ZW50cyBMUzFzYjJjdFptbHNaVDFPVDE5VFZGSWdMUzFqYkdWaGNpMXNiMmR6UFVaaGJITmxJQzB0WVc1aGJIbHphWE05YUhSMFkzQXVZMjl1Wm1sbkxtRnVZV3g1YzJselgyaDBkR053TG1GdVlXeDVjMmx6WDJoMGRHTndJQzB0ZG1WeWMybHZiajFtY21GdFpYZHZjbXRmYzNsdVl5QXRMV052Ym1acFp6MXlkVzR6WHpJd01qSmZjSEpsUlVWZmJYVjBZWFVnTFMxemFHbG1kRDF1YjIxcGJtRnNJQzB0WkdGMFlYTmxkRDFrZVY5c1pYQmZiV0ZrWjNKaGNHZ2dMUzFqWVd4cFluSmhkRzl5Y3oxdFlXbHVJQzB0YzJWc1pXTjBiM0k5YldGcGJpQXRMV05vWldOckxXWnBibWwwWlMxdmRYUndkWFE5Um1Gc2MyVWdMUzFqYUdWamF5MXZkbVZ5YkdGd2NHbHVaeTFwYm5CMWRITTlSbUZzYzJVZ0xTMWpaaTVEWVd4cFluSmhkR1ZGZG1WdWRITXRkMjl5YTJac2IzYzlhSFJqYjI1a2IzSWdMUzFqWmk1RFlXeHBZbkpoZEdWRmRtVnVkSE10ZG1WeWMybHZiajFtY21GdFpYZHZjbXRmYzNsdVl5QXRMV05tTGxKbFpIVmpaVVYyWlc1MGN5MTNiM0pyWm14dmR6MW9kR052Ym1SdmNpQXRMV05tTGxKbFpIVmpaVVYyWlc1MGN5MTJaWEp6YVc5dVBXWnlZVzFsZDI5eWExOXplVzVqSUMwdFkyWXVUV1Z5WjJWU1pXUjFZMlZrUlhabGJuUnpMWGR2Y210bWJHOTNQV2gwWTI5dVpHOXlJQzB0WTJZdVRXVnlaMlZTWldSMVkyVmtSWFpsYm5SekxYWmxjbk5wYjI0OVpuSmhiV1YzYjNKclgzTjVibU1nTFMxalppNU5aWEpuWlZObGJHVmpkR2x2YmxOMFlYUnpMWFpsY25OcGIyNDlabkpoYldWM2IzSnJYM041Ym1NZ0xTMWpaaTVRY205MmFXUmxVbVZrZFdObFpFVjJaVzUwY3kxMlpYSnphVzl1UFdaeVlXMWxkMjl5YTE5emVXNWpJQzB0Ykc5allXd3RjMk5vWldSMWJHVnlQVlJ5ZFdVPSBOdz09IDEgbm8gTFE9PVwiIn0="
    if [ -z "${render_variables}" ]; then
        >&2 echo "empty render variables"
        return "4"
    fi

    # decode
    render_variables="$( echo "${render_variables}" | base64 --decode )"

    # check files to render
    local input_files_render=( law_job_60abdeb487.sh remote_bootstrap_11d1b65a98.sh )
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
    local job_file="law_job_60abdeb487.sh"
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
