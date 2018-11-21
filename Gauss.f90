! Fonction du calcul de la solution du système AX=B selon GAUSS
Subroutine GaussMat(Jaa, Bb, Xx, Nn)

        integer, intent(IN) :: Nn
        Real, intent(IN) :: Jaa(:,:), Bb(:)
        Real, allocatable, intent(OUT) :: Xx(:)

        Real, allocatable, Dimension(:,:) :: Aa
        integer :: ii, ij, ik
        Real :: Xsum, xmult
        !character(len=40) :: A3Fmt

        ! Formation de la matrice a = [Jaa | I]
        allocate(Aa(Nn,Nn+1), Xx(Nn))

        Xx=0.
        Aa=0.
        do ii=1, Nn
            do ij=1, Nn+1
               Aa(ii,ij) = Jaa(ii,ij)
            end do
            Aa(ii, Nn+1) = Bb(ii)
        end do

        ! Elimination de Gauss
        do ik = 1, Nn-1
           do ii = ik+1, Nn
              xmult = Aa(ii, ik) / Aa(ik, ik)
              Aa(ii, ik) = 0.
              do ij = ik+1, Nn+1
                 Aa(ii, ij) = Aa(ii, ij) - xmult * Aa(ik, ij)
              end do
           end do
        end do

        Xx(Nn) = Aa(Nn, Nn+1) / Aa(Nn, Nn)
        do ii = Nn-1, 1, -1
           Xsum = 0.
           do ij= ii+1, Nn
              Xsum = Xsum + Aa(ii, ij) * Xx(ij)
           end do
           Xx(ii) = (Aa(ii, Nn+1) - Xsum) / Aa(ii, ii)
        end do

        deallocate(Aa)
 End Subroutine
