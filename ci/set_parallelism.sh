#!/bin/bash
# Export CMAKE_BUILD_PARALLEL_LEVEL and CTEST_PARALLEL_LEVEL for the following
# workflow steps, sized to the CPUs this job may use: the cgroup CPU quota when
# one is set (containers), otherwise the host CPU count, divided by
# CI_RUNNERS_PER_HOST (runner instances sharing the host, default 1).
set -eu

if command -v nproc >/dev/null 2>&1; then
    cpus=$(nproc)
else
    cpus=$(sysctl -n hw.ncpu)
fi

if [ -r /sys/fs/cgroup/cpu.max ]; then
    read -r quota period </sys/fs/cgroup/cpu.max
    if [ "$quota" != max ]; then
        quota_cpus=$((quota / period))
        if [ "$quota_cpus" -ge 1 ] && [ "$quota_cpus" -lt "$cpus" ]; then
            cpus=$quota_cpus
        fi
    fi
fi

runners=${CI_RUNNERS_PER_HOST:-1}
jobs=$((cpus / runners))
if [ "$jobs" -lt 1 ]; then
    jobs=1
fi

echo "CPU budget: $cpus CPUs, $runners runner(s) per host -> $jobs parallel jobs"
{
    echo "CMAKE_BUILD_PARALLEL_LEVEL=$jobs"
    echo "CTEST_PARALLEL_LEVEL=$jobs"
} >>"${GITHUB_ENV:?GITHUB_ENV is not set}"
