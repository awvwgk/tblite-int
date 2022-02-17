! This file is part of tblite-int.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

program integrator
   use, intrinsic :: iso_fortran_env, only : output_unit
   use mctc_env, only : wp, error_type, fatal_error, get_argument
   use mctc_io, only : structure_type, new
   use mctc_io_symbols, only : symbol_length
   use mctc_io_convert, only : aatoau
   use stdlib_io_npy, only : save_npy
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_type, only : get_cutoff
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_type, only : integral_type, new_integral
   use tblite_param, only : param_record
   use tblite_xtb_calculator, only : xtb_calculator, new_xtb_calculator
   use tblite_xtb_h0, only : get_selfenergy, get_hamiltonian
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_ipea1, only : new_ipea1_calculator
   implicit none
   character(len=*), parameter :: prog_name = "tblite-int"

   type :: int_config
      character(len=:), allocatable :: output_file
      character(len=:), allocatable :: method
      character(len=:), allocatable :: param
      real(wp), allocatable :: distance(:)
      character(len=:), allocatable :: species1, species2
      logical :: hamiltonian = .false.
      real(wp) :: accuracy = 1.0_wp
      real(wp) :: conv = aatoau
   end type int_config

   type(int_config) :: config
   type(error_type), allocatable :: error

   call get_arguments(config, error)
   call handle_error(error)

   call main(config, error)
   call handle_error(error)

contains

   subroutine main(config, error)
      type(int_config), intent(in) :: config
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(param_record) :: param
      type(xtb_calculator) :: calc
      character(len=:), allocatable :: method
      real(wp) :: cutoff
      real(wp), allocatable :: cn(:), selfenergy(:), lattr(:, :)
      type(integral_type) :: ints
      type(adjacency_list) :: list
      character(len=symbol_length) :: sym(2)
      real(wp) :: xyz(3, 2)

      sym = [character(len=symbol_length):: config%species1, config%species2]
      xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, config%distance*config%conv], [3, 2])
      call new(mol, sym, xyz)

      if (allocated(config%param)) then
         call param%load(config%param, error)
         if (.not. allocated(error)) then
            call new_xtb_calculator(calc, mol, param, error)
         end if
      else
         method = "gfn2"
         if (allocated(config%method)) method = config%method
         select case(method)
         case default
            call fatal_error(error, "Unknown method '"//method//"' requested")
         case("gfn2")
            call new_gfn2_calculator(calc, mol)
         case("gfn1")
            call new_gfn1_calculator(calc, mol)
         case("ipea1")
            call new_ipea1_calculator(calc, mol)
         end select
      end if
      if (allocated(error)) return

      if (allocated(calc%ncoord)) then
         allocate(cn(mol%nat))
         call calc%ncoord%get_cn(mol, cn)
      end if

      allocate(selfenergy(calc%bas%nsh))
      call get_selfenergy(calc%h0, mol%id, calc%bas%ish_at, calc%bas%nsh_id, cn=cn, &
         & selfenergy=selfenergy)

      cutoff = get_cutoff(calc%bas, config%accuracy)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call new_adjacency_list(list, mol, lattr, cutoff)

      call new_integral(ints, calc%bas%nao)
      call get_hamiltonian(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
         & ints%overlap, ints%dipole, ints%quadrupole, ints%hamiltonian)

      if (config%hamiltonian) then
         call save_npy(config%output_file, ints%hamiltonian)
         write(output_unit, '(a)') "[Info] Hamiltonian written to '"//config%output_file//"'"
      else
         call save_npy(config%output_file, ints%overlap)
         write(output_unit, '(a)') "[Info] Overlap written to '"//config%output_file//"'"
      end if

   end subroutine main

   subroutine get_arguments(config, error)
      type(int_config), intent(out) :: config
      type(error_type), allocatable, intent(out) :: error

      integer :: iarg, narg
      character(len=:), allocatable :: arg

      iarg = 0
      narg = command_argument_count()
      do while(iarg < narg)
         iarg = iarg + 1
         call get_argument(iarg, arg)
         select case(arg)
         case default
            if (.not.allocated(config%species1)) then
               call move_alloc(arg, config%species1)
               cycle
            end if
            if (.not.allocated(config%species2)) then
               call move_alloc(arg, config%species2)
               cycle
            end if
            if (.not.allocated(config%distance)) then
               allocate(config%distance(3))
               call get_argument_as_realv(arg, config%distance, error)
               if (allocated(error)) exit
               cycle
            end if
            call fatal_error(error, "Too many positional arguments present")
            exit

         case("--help")
            call help(output_unit)
            stop

         case("--method")
            if (allocated(config%param)) then
               call fatal_error(error, "Cannot specify method if parameter file is provided")
               exit
            end if
            iarg = iarg + 1
            call get_argument(iarg, config%method)
            if (.not.allocated(config%method)) then
               call fatal_error(error, "Missing argument for method")
               exit
            end if

         case("--param")
            if (allocated(config%param)) then
               call fatal_error(error, "Cannot specify parameter file if method is provided")
               exit
            end if
            iarg = iarg + 1
            call get_argument(iarg, config%param)
            if (.not.allocated(config%param)) then
               call fatal_error(error, "Missing argument for param")
               exit
            end if

         case("--bohr")
            config%conv = 1.0_wp

         case("-o", "--output")
            iarg = iarg + 1
            call get_argument(iarg, arg)
            if (.not.allocated(arg)) then
               call fatal_error(error, "Output file name not specified")
               exit
            end if
            call move_alloc(arg, config%output_file)

         end select
      end do

      if (.not.allocated(config%output_file)) then
         config%output_file = "tblite.npy"
      end if

      if (.not.all([allocated(config%species1), allocated(config%species2), &
         & allocated(config%distance)])) then
         call help(output_unit)
         call fatal_error(error, "Error: Missing positional arguments")
         call handle_error(error)
      end if
   end subroutine get_arguments

   subroutine get_argument_as_realv(arg, val, error)
      !> Index of command line argument, range [0:command_argument_count()]
      character(len=:), intent(in), allocatable :: arg
      !> Real value
      real(wp), intent(out) :: val(:)
      !> Error handling
      type(error_type), allocatable :: error

      integer :: stat
      integer :: i
      character(len=*), parameter :: sep = ","
      character(len=:), allocatable :: targ

      if (.not.allocated(arg)) then
         call fatal_error(error, "Cannot read real value, argument missing")
         return
      end if
      allocate(character(len=len(arg)) :: targ)
      do i = 1, len(arg)
         if (arg(i:i) == sep) then
            targ(i:i) = " "
         else
            targ(i:i) = arg(i:i)
         end if
      end do
      read(targ, *, iostat=stat) val
      if (stat /= 0) then
         call fatal_error(error, "Cannot read real value from '"//arg//"'")
         return
      end if

   end subroutine get_argument_as_realv

   subroutine help(io)
      integer, intent(in) :: io
      character(len=*), parameter :: nl = new_line('a')
      character(len=*), parameter :: help_text = &
         "Usage: "//prog_name//" [options] <element1> <element2> <distance>"//nl//&
         ""//nl//&
         "Takes two element symbols and a distance to calculate the overlap."//nl//&
         "Distance vector is given as three comma-separated values in Ångström."//nl//&
         ""//nl//&
         "Options"//nl//&
         ""//nl//&
         "      --method <name>     Parametrization of the xTB Hamiltonian to use"//nl//&
         "                          Available methods: gfn1, gfn2, ipea1 (Default: gfn2)"//nl//&
         "      --param <file>      Parametrization file to use for calculation"//nl//&
         "      --hamiltonian       Calculate Hamiltonian instead of overlap integrals"//nl//&
         "      --bohr              Use Bohr instead of Angstrom as distance unit"//nl//&
         "  -o, --output <file>     Output file for writing Hamiltonian"//nl//&
         "      --help              Show this help message"//nl//&
         ""

         write(io, '(a)') help_text
      end subroutine help

   subroutine handle_error(error)
      use, intrinsic :: iso_fortran_env, only : error_unit
      type(error_type), allocatable, intent(inout) :: error

      interface
         subroutine sys_exit(stat) bind(c, name="exit")
            use, intrinsic :: iso_c_binding, only : c_int
            integer(c_int), value :: stat
         end subroutine sys_exit
      end interface

      if (allocated(error)) then
         write(error_unit, '(a)') error%message

         ! Use exit(3) to get consistent behaviour across compilers
         call sys_exit(1)
      end if
   end subroutine handle_error

end program integrator
