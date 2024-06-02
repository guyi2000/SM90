! 数据精度定义
! ***************
module NumKind
! ***************
    implicit none
    integer(kind(1)), parameter :: ikind=kind(1), rkind=kind(0.D0)
    real(rkind), parameter :: Zero =0._rkind, One  = 1._rkind, &
            Two   = 2._rkind, Three=3._rkind, Four = 4._rkind, &
            Five  = 5._rkind, Six  =6._rkind, Seven= 7._rkind, &
            Eight = 8._rkind, Nine =9._rkind, Ten  =10._rkind, &
            Twelve=12._rkind
end module NumKind

! 结构类型定义
! ***************
module TypeDef
! ***************
    use NumKind
    implicit none

    ! 结点类型
    type :: typ_Joint
        real(rkind)     :: x, y
        integer(ikind)  :: GDOF(3)
    end type typ_Joint

    ! 单元类型
    type :: typ_Element
        integer(ikind)  :: JointNo(2), GlbDOF(6)
        real(rkind)     :: Len, CosA, SinA, EI, EA, mass
    end type typ_Element

    ! 结点荷载
    type :: typ_JointLoad
        integer(ikind)  :: JointNo, LodDOF
        real(rkind)     :: LodVal
    end type typ_JointLoad

    ! 单元荷载
    type :: typ_ElemLoad
        integer(ikind)  :: ElemNo, Indx
        real(rkind)     :: Pos, LodVal
    end type typ_ElemLoad

    contains

    ! =================================
    subroutine SetElemProp(Elem, Joint)
    ! =================================
        type(typ_Element), intent(in out) :: Elem(:)
        type(typ_Joint), intent(in)       :: Joint(:)
        
        integer(ikind) :: ie, NElem, J1, J2
        real(rkind)    :: x1, y1, x2, y2

        NElem = size(Elem)
        do ie = 1, NElem
            J1 = Elem(ie) % JointNo(1)
            J2 = Elem(ie) % JointNo(2)
            x1 = Joint(J1) % x
            y1 = Joint(J1) % y
            x2 = Joint(J2) % x
            y2 = Joint(J2) % y
            Elem(ie) % Len  = sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
            Elem(ie) % CosA = (x2 - x1) / Elem(ie) % Len
            Elem(ie) % SinA = (y2 - y1) / Elem(ie) % Len
            Elem(ie) % GlbDOF(1:3) = Joint(J1) % GDOF(:)
            Elem(ie) % GlbDOF(4:6) = Joint(J2) % GDOF(:)
        end do
        return
    end subroutine SetElemProp

    ! ====================================
    subroutine TransMatrix(ET, CosA, SinA)
    ! ====================================
        real(rkind), intent(out) :: ET(:, :)
        real(rkind), intent(in)  :: CosA, SinA
        ET = Zero
        ET(1, 1:2) = (/ CosA, SinA/)
        ET(2, 1:2) = (/-SinA, CosA/)
        if (size(ET, 1) > 2) ET(3, 3) = One
        if (size(ET, 1) > 3) ET(4:6, 4:6) = ET(1:3, 1:3)
        return
    end subroutine TransMatrix

end module TypeDef

! 变带宽模块
! ************
module BandMat
! ************
    use NumKind
    use TypeDef, only : typ_Element
    implicit none

    private
    public :: SetMatBand, DelMatBand, VarBandSolv, DiagSolv

    type, public :: typ_Kcol
        real(rkind), pointer :: row(:)
    end type typ_Kcol

    contains

    ! ===============================
    subroutine SetMatBand(Kcol, Elem)
    ! ===============================
        type(typ_Kcol), intent(in out) :: Kcol(:)
        type(typ_Element), intent(in)  :: Elem(:)
        integer(ikind)                 :: minDOF, ELocVec(6)
        integer(ikind)                 :: ie, j, NGlbDOF, NElem
        integer(ikind)                 :: row1(size(Kcol, dim=1))

        NGlbDOF = size(Kcol, 1)
        NElem   = size(Elem, 1)
        row1    = NGlbDOF

        do ie = 1, NElem
            ELocVec = Elem(ie) % GlbDOF
            minDOF  = minval(ELocVec, mask=ELocVec>0)
            where(ELocVec > 0)
                row1(ELocVec) = min(row1(ELocVec), minDOF)
            end where
        end do

        do j = 1, NGlbDOF
            allocate(Kcol(j) % row(row1(j):j))
            Kcol(j) % row = Zero
        end do
        return
    end subroutine SetMatBand

    ! ===============================
    subroutine DelMatBand(Kcol)
    ! ===============================
        type(typ_Kcol), intent(in out) :: Kcol(:)
        integer(ikind)                 :: j, NGlbDOF
        NGlbDOF = size(Kcol, 1)
        do j = 1, NGlbDOF
            deallocate(Kcol(j) % row)
        end do
        return
    end subroutine DelMatBand

    ! ===============================
    subroutine VarBandSolv(Disp, Kcol, GLoad)
    ! ===============================
        type(typ_KCol), intent(in out) :: Kcol(:)
        real(rkind),    intent(out)    :: Disp(:)
        real(rkind),    intent(in)     :: GLoad(:)
        integer(ikind)                 :: i, j, row1j, row_1, NCol
        real(rkind)                    :: s, Diag(size(Kcol, dim=1))

        NCol = size(Kcol, 1)
        Diag(:) = (/(Kcol(j) % row(j), j = 1, NCol)/)
        do j = 2, NCol
            row1j = lbound(Kcol(j) % row, 1)
            do i = row1j, j - 1
                row_1 = max(row1j, lbound(Kcol(i) % row, 1))
                s = sum(Diag(row_1:i-1) * Kcol(i) % row(row_1:i-1) &
                                        * Kcol(j) % row(row_1:i-1))
                Kcol(j)%row(i) = (Kcol(j) % row(i) - s) / Diag(i)
            end do
            s = sum(Diag(row1j:j-1) * Kcol(j) % row(row1j:j-1) ** 2)
            Diag(j) = Diag(j) - s
        end do

        Disp(:) = GLoad(:)
        do j = 2, NCol
            row1j = lbound(Kcol(j) % row, 1)
            Disp(j) = Disp(j) - sum(Kcol(j) % row(row1j:j-1) * Disp(row1j:j-1))
        end do

        Disp(:) = Disp(:) / Diag(:)
        do j = NCol, 1, -1
            row1j = lbound(Kcol(j) % row, 1)
            Disp(row1j:j-1) = Disp(row1j:j-1) - Disp(j) * Kcol(j) % row(row1j:j-1)
        end do
        return
    end subroutine VarBandSolv

    ! ===============================
    subroutine DiagSolv(Kcol, Diag)
    ! ===============================
        type(typ_Kcol), intent(in out) :: Kcol(:)
        real(rkind), intent(out)       :: Diag(:)
        integer(ikind)                 :: i, j, row1j, row_1, NCol
        real(rkind)                    :: s

        NCol = size(Kcol, 1)
        Diag(:) = (/(Kcol(j) % row(j), j = 1, NCol)/)
        do j = 2, NCol
            row1j = lbound(Kcol(j) % row, 1)
            do i = row1j, j - 1
                row_1 = max(row1j, lbound(Kcol(i) % row, 1))
                s = sum(Diag(row_1:i-1) * Kcol(i) % row(row_1:i-1) &
                                        * Kcol(j) % row(row_1:i-1))
                Kcol(j)%row(i) = (Kcol(j) % row(i) - s) / Diag(i)
            end do
            s = sum(Diag(row1j:j-1) * Kcol(j) % row(row1j:j-1) ** 2)
            Diag(j) = Diag(j) - s
        end do

        return
    end subroutine DiagSolv
end module BandMat

! 超静定结构位移求解器
! ***************
module DispMethod
! ***************
    use NumKind
    use TypeDef
    use BandMat
    implicit none

    private
    public :: SolveDisp, ElemDisp, ElemForce

    contains

    ! ===============================
    subroutine SolveDisp(Disp, Elem, Joint, JLoad, ELoad)
    ! ===============================
        real(rkind),         intent(out) :: Disp(:)
        type(typ_Element),   intent(in)  :: Elem(:)
        type(typ_Joint),     intent(in)  :: Joint(:)
        type(typ_JointLoad), intent(in)  :: JLoad(:)
        type(typ_ElemLoad),  intent(in)  :: ELoad(:)

        integer(ikind)                   :: NGlbDOF
        type(typ_Kcol), allocatable      :: Kcol(:)
        real(rkind),    allocatable      :: GLoad(:)

        Disp = Zero
        NGlbDOF = size(Disp, 1)
        allocate(GLoad(NGlbDOF))
        allocate(Kcol(NGlbDOF))

        call SetMatBand(Kcol, Elem)
        call GLoadVec(GLoad, Elem, JLoad, ELoad, Joint)
        call GStifMat(Kcol, Elem)
        call VarBandSolv(Disp, Kcol, GLoad)
        call DelMatBand(Kcol)

        return

    end subroutine SolveDisp
    
    ! ===============================
    subroutine GStifMat(Kcol, Elem)
    ! ===============================
        type(typ_Kcol),    intent(in out) :: Kcol(:)
        type(typ_Element), intent(in)     :: Elem(:)
        integer(ikind)                    :: ie, j, JGDOF, NElem
        integer(ikind)                    :: ELocVec(6)
        real(rkind)                       :: EK(6, 6), ET(6, 6)

        NElem = size(Elem, 1)
        do ie = 1, NElem
            call EStifMat(EK, Elem(ie) % Len, Elem(ie) % EI, Elem(ie) % EA)
            call TransMatrix(ET, Elem(ie) % CosA, Elem(ie) % SinA)
            EK = matmul(transpose(ET), matmul(EK, ET))
            ELocVec = Elem(ie) % GlbDOF
            do j = 1, 6
                JGDOF = ELocVec(j)
                if (JGDOF == 0) cycle
                where (ELocVec <= JGDOF .and. ELocVec > 0) 
                    Kcol(JGDOF) % row(ELocVec) = Kcol(JGDOF) % row(ELocVec) + EK(:, j)
                end where
            end do
        end do
        return
    end subroutine GStifMat

    ! ===============================
    subroutine GLoadVec(GLoad, Elem, JLoad, ELoad, Joint)
    ! ===============================
        real(rkind),         intent(out) :: GLoad(:)
        type(typ_Element),   intent(in)  :: Elem(:)
        type(typ_Joint),     intent(in)  :: Joint(:)
        type(typ_JointLoad), intent(in)  :: JLoad(:)
        type(typ_ElemLoad),  intent(in)  :: ELoad(:)

        integer(ikind)                   :: ie, NJLoad, NELoad
        real(rkind)                      :: F0(6), ET(6, 6)

        NJLoad = size(JLoad, dim=1)
        NELoad = size(ELoad, dim=1)
        GLoad(:) = Zero

        do ie = 1, NJLoad
            GLoad(Joint(JLoad(ie) % JointNo) % GDOF(JLoad(ie) % LodDOF)) = &
            GLoad(Joint(JLoad(ie) % JointNo) % GDOF(JLoad(ie) % LodDOF)) + JLoad(ie) % LodVal
        end do

        do ie = 1, NELoad
            call EFixendF(F0, ELoad(ie) % Indx, ELoad(ie) % Pos, ELoad(ie) % LodVal, Elem(ELoad(ie) % ElemNo))
            call TransMatrix(ET, Elem(ELoad(ie) % ElemNo) % CosA, Elem(ELoad(ie) % ElemNo) % SinA)
            where (Elem(ELoad(ie) % ElemNo) % GlbDOF > 0)
                GLoad(Elem(ELoad(ie) % ElemNo) % GlbDOF) = &
                GLoad(Elem(ELoad(ie) % ElemNo) % GlbDOF) - matmul(transpose(ET), F0(:))
            end where
        end do

        return
    end subroutine GLoadVec

    ! ===============================
    subroutine EStifMat(EK, ELen, EI, EA)
    ! ===============================
        real(rkind), intent(out) :: EK(6, 6)
        real(rkind), intent(in)  :: ELen, EI, EA

        real(rkind) :: EAL, EIL, EIL2, EIL3
        EAL  = EA / ELen
        EIL  = EI / ELen
        EIL2 = EI / ELen ** 2
        EIL3 = EI / ELen ** 3

        EK(1, :) = (/ EAL,           Zero,        Zero, -EAL,           Zero,        Zero/)
        EK(2, :) = (/Zero,  Twelve * EIL3,  Six * EIL2, Zero, -Twelve * EIL3,  Six * EIL2/)
        EK(3, :) = (/Zero,     Six * EIL2, Four * EIL , Zero,    -Six * EIL2,  Two * EIL /)
        EK(4, :) = (/-EAL,           Zero,        Zero,  EAL,           Zero,        Zero/)
        EK(5, :) = (/Zero, -Twelve * EIL3, -Six * EIL2, Zero,  Twelve * EIL3, -Six * EIL2/)
        EK(6, :) = (/Zero,     Six * EIL2,  Two * EIL , Zero,    -Six * EIL2, Four * EIL /)

        return
    end subroutine EStifMat

    ! ===============================
    subroutine EFixendF(F0, Indx, a, q, Elem)
    ! ===============================
        real(rkind),       intent(out) :: F0(6)
        real(rkind),       intent(in)  :: a, q
        integer(ikind),    intent(in)  :: Indx
        type(typ_Element), intent(in)  :: Elem

        real(rkind)                    :: l
        l = Elem % Len

        select case (Indx)
        case(1)
            F0(1:6)=(/Zero, -q * (a * l) / Two * (Two - Two * a ** 2 + a ** 3),              &
                            -q * (a * l) ** 2 / Twelve * (Six - Eight * a + Three * a ** 2), &
                      Zero, -q * (a * l) * a ** 2 / Two * (Two - a),                         &
                             q * (a * l) ** 2 * a / Twelve * (Four - Three * a)/)
        case(2)
            F0(1:6)=(/Zero, -q * (One - Two * a + a ** 2) * (One + Two * a), &
                            -q * (a * l) * (One - Two * a + a ** 2),         &
                      Zero, -q * a ** 2 * (Three - Two * a),                 &
                             q * a ** 2 * (One - a) * l/)
        case(3)
            if(a < One / Two) then
                F0(1) = Elem % EA * q / l
                F0(4) = -F0(1)
             else
                F0(1) = -Elem % EA * q / l
                F0(4) = -F0(1)
             end if
             F0(2) = Zero
             F0(3) = Zero
             F0(5) = Zero
             F0(6) = Zero 
        case(4)
            if(a < One / Two)then
                F0(2) =  Twelve * Elem % EI * q / l ** 3
                F0(3) =     Six * Elem % EI * q / l ** 2
                F0(5) = -F0(2)
                F0(6) =  F0(3)
             else
                F0(2) = -Twelve * Elem % EI * q / l ** 3
                F0(3) =    -Six * Elem % EI * q / l ** 2
                F0(5) = -F0(2)
                F0(6) =  F0(3)
             end if
             F0(1) = Zero
             F0(4) = Zero 
        end select

        return
    end subroutine EFixendF

    ! ===============================
    subroutine ElemDisp(EDisp, ie, Disp, Elem, ELoad)
    ! ===============================
        real(rkind),        intent(out) :: EDisp(:)
        integer(ikind),     intent(in)  :: ie
        real(rkind),        intent(in)  :: Disp(:)
        type(typ_Element),  intent(in)  :: Elem(:)
        type(typ_ElemLoad), intent(in)  :: ELoad(:)

        real(rkind)                     :: ET(6, 6)
        real(rkind)                     :: EDispLoad(6)
        integer(ikind)                  :: i
        
        do i = 1, 6
            if(Elem(ie) % GlbDOF(i) == 0) then
                EDisp(i) = Zero
            else
                EDisp(i) = Disp(Elem(ie) % GlbDOF(i))
            end if
        end do
        
        EDispLoad = Zero
        call TransMatrix(ET, Elem(ie) % CosA, Elem(ie) % SinA)
        do i = 1, size(ELoad)
            if(ELoad(i) % ElemNo == ie) then
                if(ELoad(i) % Indx == 3) then
                    if (ELoad(i) % Pos < One / Two) then
                        EDispLoad(1) = EDispLoad(1) + ELoad(i) % LodVal
                    else
                        EDispLoad(4) = EDispLoad(4) + ELoad(i) % LodVal
                    end if
                else if(ELoad(i) % Indx == 4) then
                    if (ELoad(i) % Pos < One / Two) then
                        EDispLoad(2) = EDispLoad(2) + ELoad(i) % LodVal
                    else
                        EDispLoad(5) = EDispLoad(5) + ELoad(i) % LodVal
                    end if
                end if
            end if
        end do
        EDispLoad = matmul(transpose(ET), EDispLoad)
        EDisp = EDisp + EDispLoad

        return
    end subroutine ElemDisp

    ! ===============================
    subroutine ElemForce(EForce, ie, Disp, Elem, ELoad)
    ! ===============================
        real(rkind),        intent(out) :: EForce(:)
        integer(ikind),     intent(in)  :: ie
        real(rkind),        intent(in)  :: Disp(:)
        type(typ_Element),  intent(in)  :: Elem(:)
        type(typ_ElemLoad), intent(in)  :: ELoad(:)

        real(rkind)                     :: EK(6, 6), ET(6, 6)
        real(rkind)                     :: EP(6), F0(6), EDisp(6)
        integer(ikind)                  :: NELoad, i
        NELoad = size(ELoad, 1)
        EP(:) = Zero
        do i = 1, NELoad
            if (ELoad(i) % ElemNo == ie) then
                call EFixendF(F0, ELoad(i) % Indx, ELoad(i) % Pos, &
                              ELoad(i) % LodVal, Elem(ELoad(i) % ElemNo))
                EP(:) = EP(:) - F0(:)
            end if
        end do
        call EStifMat(EK, Elem(ie) % Len, Elem(ie) % EI, Elem(ie) % EA)
        call TransMatrix(ET, Elem(ie) % CosA, Elem(ie) % SinA)

        do i = 1, 6
            if(Elem(ie) % GlbDOF(i) == 0) then
                EDisp(i) = Zero
            else
                EDisp(i) = Disp(Elem(ie) % GlbDOF(i))
            end if
        end do

        EForce(:) = matmul(EK, matmul(ET, EDisp)) - EP(:)
        EForce(1) = -EForce(1)
        return
    end subroutine ElemForce

end module DispMethod

! 自由振动频率及振型求解器
! ***************
module FreqMethod
! ***************
    use NumKind
    use TypeDef
    use BandMat

    implicit none

    private
    public :: CalculateFreq

    contains

    ! 单元动力刚度矩阵
    ! ===============================
    subroutine EDStifMat(EK, Elem, freq)
    ! ===============================
        real(rkind),       intent(out) :: EK(6, 6)
        type(typ_Element), intent(in)  :: Elem
        real(rkind),       intent(in)  :: freq
        
        real(rkind) :: ELen, EA, EI, mass
        real(rkind) :: lambda, nu, sl, cl, phi, esh, ech, inve
        real(rkind) :: B1, B2, T, R, Q, H, S, C
        real(rkind) :: EAL, EIL, EIL2, EIL3
        ELen = Elem % Len
        EA   = Elem % EA
        EI   = Elem % EI
        mass = Elem % mass

        EAL  = EA / ELen
        EIL  = EI / ELen
        EIL2 = EI / ELen ** 2
        EIL3 = EI / ELen ** 3

        nu = freq * ELen * sqrt(mass / EA)
        lambda = ELen * ((freq ** 2 * mass / EI) ** 0.25)
        
        sl   = sin(lambda)
        cl   = cos(lambda)
        inve = exp(-lambda)
        esh  = (One - inve ** 2) / Two
        ech  = (One + inve ** 2) / Two
        
        phi = inve - ech * cl
        B1  = nu / tan(nu)
        B2  = nu / sin(nu)
        T   = lambda ** 3 * (sl * ech + cl * esh) / phi
        R   = lambda ** 3 * (esh + inve * sl)     / phi
        Q   = lambda ** 2 * (esh * sl)            / phi
        H   = lambda ** 2 * (ech - inve * cl)     / phi
        S   = lambda      * (sl * ech - cl * esh) / phi
        C   = lambda      * (esh - inve * sl)     / phi

        EK(1, :) = (/ B1 * EAL,      Zero,      Zero, -B2 * EAL,      Zero,      Zero/)
        EK(2, :) = (/     Zero,  T * EIL3,  Q * EIL2,      Zero, -R * EIL3,  H * EIL2/)
        EK(3, :) = (/     Zero,  Q * EIL2,  S * EIL ,      Zero, -H * EIL2,  C * EIL /)
        EK(4, :) = (/-B2 * EAL,      Zero,      Zero,  B1 * EAL,      Zero,      Zero/)
        EK(5, :) = (/     Zero, -R * EIL3, -H * EIL2,      Zero,  T * EIL3, -Q * EIL2/)
        EK(6, :) = (/     Zero,  H * EIL2,  C * EIL ,      Zero, -Q * EIL2,  S * EIL /)

        return
    end subroutine EDStifMat

    ! ===============================
    subroutine GDStifMat(Kcol, Elem, freq)
    ! ===============================
        type(typ_Kcol),    intent(in out) :: Kcol(:)
        type(typ_Element), intent(in)     :: Elem(:)
        real(rkind),       intent(in)     :: freq

        integer(ikind)                    :: ie, j, JGDOF, NElem
        integer(ikind)                    :: ELocVec(6)
        real(rkind)                       :: EK(6, 6), ET(6, 6)

        NElem = size(Elem, 1)
        do ie = 1, NElem
            call EDStifMat(EK, Elem(ie), freq)
            call TransMatrix(ET, Elem(ie) % CosA, Elem(ie) % SinA)
            EK = matmul(transpose(ET), matmul(EK, ET))
            ELocVec = Elem(ie) % GlbDOF
            do j = 1, 6
                JGDOF = ELocVec(j)
                if (JGDOF == 0) cycle
                where (ELocVec <= JGDOF .and. ELocVec > 0) 
                    Kcol(JGDOF) % row(ELocVec) = Kcol(JGDOF) % row(ELocVec) + EK(:, j)
                end where
            end do
        end do
        return
    end subroutine GDStifMat

    ! ===============================
    subroutine Calculate_J0(J0, freq, Elem)
    ! ===============================
        integer(ikind),    intent(out) :: J0
        real(rkind),       intent(in)  :: freq
        type(typ_Element), intent(in)  :: Elem(:)

        integer(ikind)                 :: NElem, Ja, Jb, ie, sg
        real(rkind)                    :: nu, lambda, inve

        NElem = size(Elem)
        J0 = 0
        do ie = 1, NElem
            nu = freq * Elem(ie) % Len * sqrt(Elem(ie) % mass / Elem(ie) % EA)
            lambda = Elem(ie) % Len * ((freq ** 2 * Elem(ie) % mass / Elem(ie) % EI) ** 0.25)

            Ja   = int(nu / acos(-One))
            inve = exp(-lambda)
            sg   = int(sign(One, inve - cos(lambda) * (One + inve ** 2) / Two))
            Jb   = int(lambda / acos(-One)) - (1 - (-1) ** int(lambda / acos(-One)) * sg) / 2

            J0 = J0 + Ja + Jb
        end do
        return
    end subroutine Calculate_J0

    ! ===============================
    subroutine Calculate_JK(Jk, freq, Kcol, Elem)
    ! ===============================
        integer(ikind),    intent(out)    :: Jk
        real(rkind),       intent(in)     :: freq
        type(typ_Kcol),    intent(in out) :: Kcol(:)
        type(typ_Element), intent(in)     :: Elem(:)

        real(rkind), allocatable          :: Diag(:)
        integer(ikind)                    :: NGlbDOF
        
        NGlbDOF = size(Kcol)
        allocate(Diag(NGlbDOF))

        call SetMatBand(Kcol, Elem)
        call GDStifMat(Kcol, Elem, freq)
        call DiagSolv(Kcol, Diag)
        Jk = count(mask=Diag < Zero, dim=1)

        call DelMatBand(Kcol)
        deallocate(Diag)
        return
    end subroutine Calculate_JK

    ! ===============================
    subroutine Calculate_kFreq(freq, kfreq, Toler, GKcol, Elem)
    ! ===============================
        integer(ikind),    intent(in)     :: kfreq
        real(rkind),       intent(in)     :: Toler
        type(typ_Kcol),    intent(in out) :: GKcol(:)
        type(typ_Element), intent(in)     :: Elem(:)
        real(rkind),       intent(out)    :: freq

        integer(ikind)                    :: J0, Jk
        real(rkind)                       :: freq1, freq2

        freq1 = One
        freq2 = Ten
        do 
            call Calculate_J0(J0, freq1, Elem)
            call Calculate_JK(Jk, freq1, GKcol, Elem)
            J0 = J0 + Jk
            if (J0 < kfreq) exit
            freq1 = freq1 / Two
        end do
        do 
            call Calculate_J0(J0, freq2, Elem)
            call Calculate_Jk(Jk, freq2, GKcol, Elem)
            J0 = J0 + Jk
            if (J0 > kfreq) exit
            freq2 = freq2 * Two
        end do

        do
            freq = (freq1 + freq2) / Two
            call Calculate_J0(J0, freq, Elem)
            call Calculate_JK(Jk, freq, GKcol, Elem)
            J0 = J0 + Jk
            if (J0 >= kfreq) then
                freq2 = freq
            else 
                freq1 = freq
            end if
            if((freq2 - freq1) <= Toler * (One + freq2)) exit
        end do
        freq = (freq1 + freq2) / Two
        return
    end subroutine Calculate_kFreq

    ! ===============================
    subroutine CalculateFreq(Freq, FreqStart, Toler, Elem, NGlbDOF)
    ! ===============================
        integer(ikind),    intent(in)  :: FreqStart, NGlbDOF
        real(rkind),       intent(in)  :: Toler
        type(typ_Element), intent(in)  :: Elem(:)
        real(rkind),       intent(out) :: Freq(:)

        type(typ_Kcol), allocatable    :: GKcol(:)
        integer(ikind)                 :: k, NFreq

        NFreq = size(Freq)
        allocate(GKcol(NGlbDOF))

        do k = FreqStart, FreqStart + NFreq - 1
            call Calculate_kFreq(freq(k), k, Toler, GKcol, Elem)
        end do

        return
    end subroutine CalculateFreq
end module FreqMethod

! 主程序
! =====================================
program SM_90
! =====================================
    use NumKind
    use TypeDef
    use DispMethod
    use FreqMethod

    implicit none

    integer(ikind)                   :: ProbType
    integer(ikind)                   :: NElem, NJoint, NGlbDOF
    integer(ikind)                   :: NJLoad, NELoad
    integer(ikind)                   :: NFreq, FreqStart
    real(rkind)                      :: Tol
    type(typ_Joint),     allocatable :: Joint(:)
    type(typ_Element),   allocatable :: Elem(:)
    type(typ_JointLoad), allocatable :: JLoad(:)
    type(typ_ElemLoad),  allocatable :: ELoad(:)
    real(rkind),         allocatable :: Disp(:)
    real(rkind),         allocatable :: Freq(:)

    integer(ikind)                   :: i, ie

    call InputData()
    call SetElemProp(Elem, Joint)
    if (ProbType == 1) then
        call SolveDisp(Disp, Elem, Joint, JLoad, ELoad)
    else if(ProbType == 2) then
        call CalculateFreq(Freq, FreqStart, Tol, Elem, NGlbDOF)
    end if
    call OutputResults()

    stop

    contains
    ! --------------------
    subroutine InputData()
    ! --------------------
        open(3, file="SM90.IPT", status="OLD", position="REWIND", &
             action="read")

        read(3, *) ProbType
        if (ProbType == 1) then
            read(3, *) NElem, NJoint, NGlbDOF, NJLoad, NELoad

            allocate(Joint(NJoint))
            allocate(Elem(NElem))
            allocate(JLoad(NJLoad))
            allocate(ELoad(NELoad))
            allocate(Disp(NGlbDOF))

            Disp = Zero

            read(3, *) (Joint(i), i = 1, NJoint)
            read(3, *) (Elem(ie) % JointNo, Elem(ie) % EA, Elem(ie) % EI, ie = 1, NElem)
            if (NJLoad > 0) read(3, *) (JLoad(i), i = 1, NJLoad)
            if (NELoad > 0) read(3, *) (ELoad(i), i = 1, NELoad)
        else if (ProbType == 2) then
            read(3, *) NFreq, FreqStart, Tol
            read(3, *) NElem, NJoint, NGlbDOF, NJLoad, NELoad

            allocate(Joint(NJoint))
            allocate(Elem(NElem))
            allocate(Freq(FreqStart:FreqStart+NFreq-1))

            read(3, *) (Joint(i), i = 1, NJoint)
            read(3, *) (Elem(ie) % JointNo, Elem(ie) % EA, Elem(ie) % EI, Elem(ie) % mass, ie = 1, NElem)
        else 
            return
        end if
        close(3)
        return
    end subroutine InputData

    ! ------------------------
    subroutine OutputResults()
    ! ------------------------
        real(rkind) :: EDisp(6), EForce(6)

        open(4, file="SMCAI90.OUT", action="write")
        if (ProbType == 1) then
            write(4, *) 10, 0
            do ie = 1, size(Elem)
                call ElemDisp(EDisp, ie, Disp, Elem, ELoad)
                write(4, *) EDisp
            end do
            do ie = 1, size(Elem)
                call ElemForce(EForce, ie, Disp, Elem, ELoad)
                write(4, *) EForce
            end do
        else if (ProbType == 2) then
            write(4, *) 10, 0
            write(4, *) NFreq
            do i = FreqStart, FreqStart + NFreq - 1
                write(4, *) Freq(i), 0
                do ie = 1, NElem
                    write(4, *) 0, 0, 0, 0, 0, 0
                end do
            end do
        else 
            return
        end if
        close(4)
        return

    end subroutine OutputResults

end program SM_90