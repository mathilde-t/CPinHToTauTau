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
        ['1']="columnflow.tasks.reduction ReduceEvents LS1sb2ctZmlsZT1OT19TVFIgLS1jbGVhci1sb2dzPUZhbHNlIC0tYW5hbHlzaXM9aHR0Y3AuY29uZmlnLmFuYWx5c2lzX2h0dGNwLmFuYWx5c2lzX2h0dGNwIC0tdmVyc2lvbj1mcmFtZXdvcmtfc3luYyAtLWNvbmZpZz1ydW4zXzIwMjJfcHJlRUVfbXV0YXVfbGltaXRlZCAtLXNoaWZ0PW5vbWluYWwgLS1kYXRhc2V0PWhfZ2dmX2h0dF9jcG9fZmlsdGVyZWQgLS1jYWxpYnJhdG9ycz1tYWluIC0tc2VsZWN0b3I9bWFpbiAtLXNlbGVjdG9yLXN0ZXBzPScnIC0tY2hlY2stZmluaXRlLW91dHB1dD1GYWxzZSAtLWNoZWNrLW92ZXJsYXBwaW5nLWlucHV0cz1GYWxzZSAtLWNmLkNhbGlicmF0ZUV2ZW50cy13b3JrZmxvdz1odGNvbmRvciAtLWNmLkNhbGlicmF0ZUV2ZW50cy12ZXJzaW9uPWZyYW1ld29ya19zeW5jIC0tY2YuU2VsZWN0RXZlbnRzLXdvcmtmbG93PWxvY2FsIC0tY2YuU2VsZWN0RXZlbnRzLXZlcnNpb249ZnJhbWV3b3JrX3N5bmMgLS1jZi5NZXJnZVJlZHVjZWRFdmVudHMtd29ya2Zsb3c9aHRjb25kb3IgLS1jZi5NZXJnZVJlZHVjZWRFdmVudHMtdmVyc2lvbj1mcmFtZXdvcmtfc3luYyAtLWNmLk1lcmdlU2VsZWN0aW9uU3RhdHMtdmVyc2lvbj1mcmFtZXdvcmtfc3luYyAtLWNmLlByb3ZpZGVSZWR1Y2VkRXZlbnRzLXZlcnNpb249ZnJhbWV3b3JrX3N5bmMgLS1sb2NhbC1zY2hlZHVsZXI9VHJ1ZQ== MA== 1 no LQ=="
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
    local render_variables="eyJsYXdfY29uZmlnX2ZpbGUiOiAiJENGX1JFUE9fQkFTRS9sYXcuY2ZnIiwgImNmX3JlcG9fdXJpcyI6ICJyb290Oi8vZW9zdXNlci5jZXJuLmNoLy9lb3MvcHJvamVjdC9pL2lwaGN0YXUvcHVibGljL213aXR0L0NQaW5IVG9UYXVUYXVPdXRwdXQvaGFtYm91cmcvY2Zfc3RvcmUvYW5hbHlzaXNfaHR0Y3AvY2YuQnVuZGxlUmVwbyIsICJjZl9yZXBvX3BhdHRlcm4iOiAiQ1BpbkhUb1RhdVRhdS5jZWM2MWI3YTdiMzI1NzFmYmUyNzg5OWJiNTRkYmM1MWFmMWJlY2M1LlteXFwuXSsudGd6IiwgImNmX3NvZnR3YXJlX3VyaXMiOiAicm9vdDovL2Vvc3VzZXIuY2Vybi5jaC8vZW9zL3Byb2plY3QvaS9pcGhjdGF1L3B1YmxpYy9td2l0dC9DUGluSFRvVGF1VGF1T3V0cHV0L2hhbWJvdXJnL2NmX3N0b3JlL2FuYWx5c2lzX2h0dGNwL2NmLkJ1bmRsZVNvZnR3YXJlIiwgImNmX3NvZnR3YXJlX3BhdHRlcm4iOiAic29mdHdhcmUuW15cXC5dKy50Z3oiLCAiY2ZfYmFzaF9zYW5kYm94X3VyaXMiOiAiXCJyb290Oi8vZW9zdXNlci5jZXJuLmNoLy9lb3MvcHJvamVjdC9pL2lwaGN0YXUvcHVibGljL213aXR0L0NQaW5IVG9UYXVUYXVPdXRwdXQvaGFtYm91cmcvY2Zfc3RvcmUvYW5hbHlzaXNfaHR0Y3AvY2YuQnVuZGxlQmFzaFNhbmRib3hcIiBcInJvb3Q6Ly9lb3N1c2VyLmNlcm4uY2gvL2Vvcy9wcm9qZWN0L2kvaXBoY3RhdS9wdWJsaWMvbXdpdHQvQ1BpbkhUb1RhdVRhdU91dHB1dC9oYW1ib3VyZy9jZl9zdG9yZS9hbmFseXNpc19odHRjcC9jZi5CdW5kbGVCYXNoU2FuZGJveFwiIFwicm9vdDovL2Vvc3VzZXIuY2Vybi5jaC8vZW9zL3Byb2plY3QvaS9pcGhjdGF1L3B1YmxpYy9td2l0dC9DUGluSFRvVGF1VGF1T3V0cHV0L2hhbWJvdXJnL2NmX3N0b3JlL2FuYWx5c2lzX2h0dGNwL2NmLkJ1bmRsZUJhc2hTYW5kYm94XCIiLCAiY2ZfYmFzaF9zYW5kYm94X3BhdHRlcm5zIjogIlwiY2ZfNWRlNzIxM2MuMWM1NGY4NjY2MS5bXlxcLl0rLnRnelwiIFwidmVudl9jb2x1bW5hcl80ZWY0NTY3NC4wN2Y3MTBhMmQ3LlteXFwuXSsudGd6XCIgXCJ2ZW52X2NvbHVtbmFyX3hnYl9mMDUxZjcyZC5lZWVkY2Y4MWI4LlteXFwuXSsudGd6XCIiLCAiY2ZfYmFzaF9zYW5kYm94X25hbWVzIjogIlwiY2ZcIiBcInZlbnZfY29sdW1uYXJcIiBcInZlbnZfY29sdW1uYXJfeGdiXCIiLCAiY2ZfYm9vdHN0cmFwX25hbWUiOiAiaHRjb25kb3Jfc3RhbmRhbG9uZSIsICJjZl9odGNvbmRvcl9mbGF2b3IiOiAiY2VybiIsICJjZl9wcmVfc2V0dXBfY29tbWFuZCI6ICIiLCAiY2ZfcG9zdF9zZXR1cF9jb21tYW5kIjogIiIsICJjZl9yZW1vdGVfbGNnX3NldHVwIjogIi9jdm1mcy9ncmlkLmNlcm4uY2gvYWxtYTktdWktdGVzdC9ldGMvcHJvZmlsZS5kL3NldHVwLWFsbWE5LXRlc3Quc2giLCAiY2ZfcmVtb3RlX2xjZ19zZXR1cF9mb3JjZSI6ICIiLCAiY2ZfYmFzZSI6ICIvYWZzL2Nlcm4uY2gvdXNlci9tL213aXR0L3B1YmxpYy9DUGluSFRvVGF1VGF1L21vZHVsZXMvY29sdW1uZmxvdyIsICJjZl9yZXBvX2Jhc2UiOiAiL2Fmcy9jZXJuLmNoL3VzZXIvbS9td2l0dC9wdWJsaWMvQ1BpbkhUb1RhdVRhdSIsICJjZl9jZXJuX3VzZXIiOiAibXdpdHQiLCAiY2Zfc3RvcmVfbmFtZSI6ICJjZl9zdG9yZSIsICJjZl9zdG9yZV9sb2NhbCI6ICIvYWZzL2Nlcm4uY2gvdXNlci9tL213aXR0L3B1YmxpYy9DUGluSFRvVGF1VGF1L2RhdGEvY2Zfc3RvcmUiLCAiY2ZfbG9jYWxfc2NoZWR1bGVyIjogInRydWUiLCAiam9iX2ZpbGUiOiAibGF3X2pvYl9kMDg1OTliYzA4LnNoIiwgImV4ZWN1dGFibGVfZmlsZSI6ICJodGNvbmRvcl93cmFwcGVyX2RhZThhNzM5MGEuc2giLCAiYm9vdHN0cmFwX2ZpbGUiOiAicmVtb3RlX2Jvb3RzdHJhcF9kMDdlNTBiYTk1LnNoIiwgInZvbXNwcm94eV9maWxlIjogIng1MDl1cF91MTUxNzkwX2E4ZjEyZmQxZjUiLCAid2xjZ190b29scyI6ICJsYXdfd2xjZ190b29sc180MDFkOWI5MDU3LnNoIiwgImlucHV0X2ZpbGVzIjogImxhd19qb2JfZDA4NTk5YmMwOC5zaCBodGNvbmRvcl93cmFwcGVyX2RhZThhNzM5MGEuc2ggcmVtb3RlX2Jvb3RzdHJhcF9kMDdlNTBiYTk1LnNoIHg1MDl1cF91MTUxNzkwX2E4ZjEyZmQxZjUgbGF3X3dsY2dfdG9vbHNfNDAxZDliOTA1Ny5zaCIsICJpbnB1dF9maWxlc19yZW5kZXIiOiAibGF3X2pvYl9kMDg1OTliYzA4LnNoIHJlbW90ZV9ib290c3RyYXBfZDA3ZTUwYmE5NS5zaCIsICJsb2dfZmlsZSI6ICJzdGRhbGwkKGxhd19qb2JfcG9zdGZpeCkudHh0IiwgImh0Y29uZG9yX2pvYl9hcmd1bWVudHNfbWFwIjogIlsnMSddPVwiY29sdW1uZmxvdy50YXNrcy5yZWR1Y3Rpb24gUmVkdWNlRXZlbnRzIExTMXNiMmN0Wm1sc1pUMU9UMTlUVkZJZ0xTMWpiR1ZoY2kxc2IyZHpQVVpoYkhObElDMHRZVzVoYkhsemFYTTlhSFIwWTNBdVkyOXVabWxuTG1GdVlXeDVjMmx6WDJoMGRHTndMbUZ1WVd4NWMybHpYMmgwZEdOd0lDMHRkbVZ5YzJsdmJqMW1jbUZ0WlhkdmNtdGZjM2x1WXlBdExXTnZibVpwWnoxeWRXNHpYekl3TWpKZmNISmxSVVZmYlhWMFlYVmZiR2x0YVhSbFpDQXRMWE5vYVdaMFBXNXZiV2x1WVd3Z0xTMWtZWFJoYzJWMFBXaGZaMmRtWDJoMGRGOWpjRzlmWm1sc2RHVnlaV1FnTFMxallXeHBZbkpoZEc5eWN6MXRZV2x1SUMwdGMyVnNaV04wYjNJOWJXRnBiaUF0TFhObGJHVmpkRzl5TFhOMFpYQnpQU2NuSUMwdFkyaGxZMnN0Wm1sdWFYUmxMVzkxZEhCMWREMUdZV3h6WlNBdExXTm9aV05yTFc5MlpYSnNZWEJ3YVc1bkxXbHVjSFYwY3oxR1lXeHpaU0F0TFdObUxrTmhiR2xpY21GMFpVVjJaVzUwY3kxM2IzSnJabXh2ZHoxb2RHTnZibVJ2Y2lBdExXTm1Ma05oYkdsaWNtRjBaVVYyWlc1MGN5MTJaWEp6YVc5dVBXWnlZVzFsZDI5eWExOXplVzVqSUMwdFkyWXVVMlZzWldOMFJYWmxiblJ6TFhkdmNtdG1iRzkzUFd4dlkyRnNJQzB0WTJZdVUyVnNaV04wUlhabGJuUnpMWFpsY25OcGIyNDlabkpoYldWM2IzSnJYM041Ym1NZ0xTMWpaaTVOWlhKblpWSmxaSFZqWldSRmRtVnVkSE10ZDI5eWEyWnNiM2M5YUhSamIyNWtiM0lnTFMxalppNU5aWEpuWlZKbFpIVmpaV1JGZG1WdWRITXRkbVZ5YzJsdmJqMW1jbUZ0WlhkdmNtdGZjM2x1WXlBdExXTm1MazFsY21kbFUyVnNaV04wYVc5dVUzUmhkSE10ZG1WeWMybHZiajFtY21GdFpYZHZjbXRmYzNsdVl5QXRMV05tTGxCeWIzWnBaR1ZTWldSMVkyVmtSWFpsYm5SekxYWmxjbk5wYjI0OVpuSmhiV1YzYjNKclgzTjVibU1nTFMxc2IyTmhiQzF6WTJobFpIVnNaWEk5VkhKMVpRPT0gTUE9PSAxIG5vIExRPT1cIiJ9"
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
