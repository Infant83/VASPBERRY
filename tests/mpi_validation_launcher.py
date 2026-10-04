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
