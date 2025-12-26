module interface
   
    interface
        subroutine zheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork, lrwork, liwork
            complex(8), intent(inout) :: a(lda, *)
            real(8), intent(out) :: w(*)
            complex(8), intent(inout) :: work(*)
            real(8), intent(inout) :: rwork(*)
            integer, intent(inout) :: iwork(*)
            integer, intent(out) :: info
        end subroutine zheevd

        subroutine zheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &
                          m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, &
                          iwork, liwork, info)
            character*1, intent(in) :: jobz, range, uplo
            integer, intent(in) :: n, lda, ldz, lwork, lrwork, liwork
            complex(8), intent(inout) :: a(lda, *)
            real(8), intent(in) :: vl, vu, abstol
            integer, intent(in) :: il, iu
            integer, intent(out) :: m, isuppz(*), iwork(*), info
            real(8), intent(out) :: w(*)
            complex(8), intent(out) :: z(ldz, *)
            real(8), intent(out) :: rwork(*)
            complex(8), intent(out) :: work(*)
        end subroutine zheevr

        subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            character*1, intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            complex(8), intent(in) :: alpha, beta
            complex(8), intent(in) :: a(lda, *), b(ldb, *)
            complex(8), intent(inout) :: c(ldc, *)
        end subroutine zgemm

        subroutine dgetrf(m, n, aa, lda, ipiv, info)
            integer, intent(in) :: m, n, lda
            integer, intent(out) :: ipiv(*)
            integer, intent(out) :: info
            real(8), intent(inout) :: aa(lda, *)
        end subroutine dgetrf

        subroutine dgetri(n, aa, lda, ipiv, work, lwork, info)
            integer, intent(in) :: n, lda, lwork
            integer, intent(in) :: ipiv(*)
            integer, intent(out) :: info
            real(8), intent(inout) :: aa(lda, *)
            real(8), intent(inout) :: work(*)
        end subroutine dgetri
    end interface


    
end module interface