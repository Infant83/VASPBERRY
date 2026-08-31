program test_z2_finalize_guard
  implicit none
  character(75) :: legacy_base

  call get_command_argument(1, legacy_base)
  if (len_trim(legacy_base) == 0) legacy_base = "NFIELD"
  call finalize_z2_pass(trim(legacy_base))
end program test_z2_finalize_guard
