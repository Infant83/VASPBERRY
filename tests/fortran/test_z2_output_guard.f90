program test_z2_output_guard
  implicit none
  character(4096) :: input_file, legacy_base

  call get_command_argument(1, input_file)
  call get_command_argument(2, legacy_base)
  if (len_trim(input_file) == 0) stop 2
  if (len_trim(legacy_base) == 0) legacy_base = "NFIELD"
  call begin_z2_field_outputs(trim(input_file), trim(legacy_base))
end program test_z2_output_guard
