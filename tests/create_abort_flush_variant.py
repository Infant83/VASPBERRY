"""Create a separate diagnostic source; leave the baseline source untouched."""
from pathlib import Path
import hashlib
import json

source = Path('vaspberry.f')
target = Path('vaspberry-abort-flush.f')
data = source.read_bytes()
old = b"      subroutine vaspberry_fail\r\n#ifdef MPI_USE\r\n      include 'mpif.h'\r\n      integer ierr\r\n      call MPI_ABORT(MPI_COMM_WORLD,1,ierr)\r\n#endif\r\n      stop 1\r\n      end subroutine vaspberry_fail\r\n"
new = b"      subroutine vaspberry_fail\r\n#ifdef MPI_USE\r\n      include 'mpif.h'\r\n      integer ierr\r\n#endif\r\n      integer ios\r\nc     Best-effort local diagnostic flush; no collective on an error path.\r\n      flush(0,iostat=ios)\r\n      flush(6,iostat=ios)\r\n#ifdef MPI_USE\r\n      call MPI_ABORT(MPI_COMM_WORLD,1,ierr)\r\n#endif\r\n      stop 1\r\n      end subroutine vaspberry_fail\r\n"
assert data.count(old) == 1
assert not target.exists()
target.write_bytes(data.replace(old, new))
print(json.dumps({'source_sha256': hashlib.sha256(data).hexdigest(),
                  'variant_sha256': hashlib.sha256(target.read_bytes()).hexdigest(),
                  'change': 'local unit-0/unit-6 flush before unchanged MPI_ABORT'}))
