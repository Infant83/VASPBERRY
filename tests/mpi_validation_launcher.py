"""Explicit launch policy for MPI validation; never translates exit status."""


def mpi_validation_command(command, *, ignore_sigpipe=False):
    """Optionally retain ignored SIGPIPE across exec for Intel Hydra aborts.

    The fixed shell script receives the original argv as separate arguments.
    exec leaves the MPI launcher's real exit/signal status visible to the
    existing validation guards. Other signals and the native binaries are
    untouched; callers opt in only for the diagnosed Intel environment.
    """
    if not ignore_sigpipe:
        return list(command)
    return ["/bin/sh", "-c", 'trap "" PIPE; exec "$@"',
            "vaspberry-mpi-validation", *command]


def mpi_payload_command(command, *, rank_logs=False):
    """Capture local MPI process streams without the launcher's I/O relay.

    This opt-in is for single-host validation. Each shell's process ID gives
    its files a unique name inside the invocation's fresh working directory.
    Fixed shell code and exec preserve literal argv and the payload exit code.
    """
    if not rank_logs:
        return list(command)
    return ["/bin/sh", "-c",
            'exec "$@" >"mpi-process-$$.stdout.log" 2>"mpi-process-$$.stderr.log"',
            "vaspberry-mpi-process-logs", *command]


def mpi_process_logs(directory, *, stream=None):
    """Return only actual process logs; optionally preserve stderr-only checks."""
    streams = ("stdout", "stderr") if stream is None else (stream,)
    if any(value not in ("stdout", "stderr") for value in streams):
        raise ValueError("invalid process log stream")
    return sorted(path for value in streams
                  for path in directory.glob(f"mpi-process-*.{value}.log")
                  if path.is_file())
