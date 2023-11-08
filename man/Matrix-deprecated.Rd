\name{Matrix-deprecated}
\title{Deprecated Functions in Package \pkg{Matrix}}
%
\keyword{internal}
\keyword{misc}
%
\alias{Matrix-deprecated}
%
\alias{..2dge}
\alias{.C2nC}
\alias{.SuiteSparse_version}
\alias{.T2Cmat}
\alias{.asmatrix}
\alias{.dense2sy}
\alias{.diag2mat}
\alias{.diag2sT}
\alias{.diag2tT}
\alias{.dsy2dsp}
\alias{.dsy2mat}
\alias{.dxC2mat}
\alias{.m2dgC}
\alias{.m2lgC}
\alias{.m2ngCn}
\alias{.m2ngC}
\alias{.m2ngTn}
\alias{.n2dgT}
\alias{.nC2d}
\alias{.nC2l}
%
\alias{chol2inv,CHMfactor-method}
\alias{coerce,CHMfactor,Matrix-method}
\alias{coerce,CHMfactor,CsparseMatrix-method}
\alias{coerce,CHMfactor,RsparseMatrix-method}
\alias{coerce,CHMfactor,TsparseMatrix-method}
\alias{coerce,CHMfactor,dMatrix-method}
\alias{coerce,CHMfactor,dsparseMatrix-method}
\alias{coerce,CHMfactor,pMatrix-method}
\alias{coerce,CHMfactor,sparseMatrix-method}
\alias{coerce,CHMfactor,triangularMatrix-method}
\alias{coerce,RsparseMatrix,dgeMatrix-method}
\alias{coerce,ddenseMatrix,dgeMatrix-method}
\alias{coerce,ddiMatrix,dgCMatrix-method}
\alias{coerce,ddiMatrix,dgeMatrix-method}
\alias{coerce,ddiMatrix,dtCMatrix-method}
\alias{coerce,dgCMatrix,dgTMatrix-method}
\alias{coerce,dgCMatrix,dgeMatrix-method}
\alias{coerce,dgCMatrix,dsCMatrix-method}
\alias{coerce,dgCMatrix,dtCMatrix-method}
\alias{coerce,dgCMatrix,lgCMatrix-method}
\alias{coerce,dgCMatrix,ngCMatrix-method}
\alias{coerce,dgTMatrix,dgCMatrix-method}
\alias{coerce,dgTMatrix,dgeMatrix-method}
\alias{coerce,dgTMatrix,dsTMatrix-method}
\alias{coerce,dgTMatrix,dtCMatrix-method}
\alias{coerce,dgTMatrix,dtTMatrix-method}
\alias{coerce,dgeMatrix,dgCMatrix-method}
\alias{coerce,dgeMatrix,dgTMatrix-method}
\alias{coerce,dgeMatrix,dsTMatrix-method}
\alias{coerce,dgeMatrix,dspMatrix-method}
\alias{coerce,dgeMatrix,dsyMatrix-method}
\alias{coerce,dgeMatrix,dtrMatrix-method}
\alias{coerce,dgeMatrix,lgeMatrix-method}
\alias{coerce,dsCMatrix,dgCMatrix-method}
\alias{coerce,dsCMatrix,dgTMatrix-method}
\alias{coerce,dsCMatrix,dgeMatrix-method}
\alias{coerce,dsCMatrix,dsRMatrix-method}
\alias{coerce,dsCMatrix,dsTMatrix-method}
\alias{coerce,dsCMatrix,dsyMatrix-method}
\alias{coerce,dsCMatrix,lsCMatrix-method}
\alias{coerce,dsCMatrix,nsCMatrix-method}
\alias{coerce,dsTMatrix,dgTMatrix-method}
\alias{coerce,dsTMatrix,dgeMatrix-method}
\alias{coerce,dsTMatrix,dsCMatrix-method}
\alias{coerce,dsTMatrix,dsyMatrix-method}
\alias{coerce,dsTMatrix,lsTMatrix-method}
\alias{coerce,dspMatrix,dsyMatrix-method}
\alias{coerce,dspMatrix,lspMatrix-method}
\alias{coerce,dsyMatrix,dsCMatrix-method}
\alias{coerce,dsyMatrix,dsTMatrix-method}
\alias{coerce,dsyMatrix,dspMatrix-method}
\alias{coerce,dsyMatrix,lsyMatrix-method}
\alias{coerce,dtCMatrix,dgCMatrix-method}
\alias{coerce,dtCMatrix,dgTMatrix-method}
\alias{coerce,dtCMatrix,dgeMatrix-method}
\alias{coerce,dtCMatrix,dsCMatrix-method}
\alias{coerce,dtCMatrix,dtTMatrix-method}
\alias{coerce,dtCMatrix,dtrMatrix-method}
\alias{coerce,dtCMatrix,ltCMatrix-method}
\alias{coerce,dtCMatrix,ntCMatrix-method}
\alias{coerce,dtTMatrix,dgTMatrix-method}
\alias{coerce,dtTMatrix,dgeMatrix-method}
\alias{coerce,dtTMatrix,dtCMatrix-method}
\alias{coerce,dtTMatrix,dtrMatrix-method}
\alias{coerce,dtpMatrix,dtTMatrix-method}
\alias{coerce,dtpMatrix,dtrMatrix-method}
\alias{coerce,dtpMatrix,ltpMatrix-method}
\alias{coerce,dtrMatrix,dtpMatrix-method}
\alias{coerce,dtrMatrix,ltrMatrix-method}
\alias{coerce,indMatrix,ngTMatrix-method}
\alias{coerce,indMatrix,ngeMatrix-method}
\alias{coerce,lMatrix,dgCMatrix-method}
\alias{coerce,lgCMatrix,dgCMatrix-method}
\alias{coerce,lgCMatrix,lgTMatrix-method}
\alias{coerce,lgCMatrix,lgeMatrix-method}
\alias{coerce,lgCMatrix,ltCMatrix-method}
\alias{coerce,lgTMatrix,dgTMatrix-method}
\alias{coerce,lgTMatrix,lgCMatrix-method}
\alias{coerce,lgTMatrix,lgeMatrix-method}
\alias{coerce,lgTMatrix,lsCMatrix-method}
\alias{coerce,lgTMatrix,ltTMatrix-method}
\alias{coerce,lgeMatrix,dgeMatrix-method}
\alias{coerce,lgeMatrix,lgCMatrix-method}
\alias{coerce,lgeMatrix,lgTMatrix-method}
\alias{coerce,lgeMatrix,lspMatrix-method}
\alias{coerce,lgeMatrix,lsyMatrix-method}
\alias{coerce,lgeMatrix,ltpMatrix-method}
\alias{coerce,lgeMatrix,ltrMatrix-method}
\alias{coerce,lsCMatrix,dsCMatrix-method}
\alias{coerce,lsCMatrix,lgCMatrix-method}
\alias{coerce,lsCMatrix,lgTMatrix-method}
\alias{coerce,lsCMatrix,lsTMatrix-method}
\alias{coerce,lsTMatrix,lgCMatrix-method}
\alias{coerce,lsTMatrix,lgTMatrix-method}
\alias{coerce,lsTMatrix,lsCMatrix-method}
\alias{coerce,lsTMatrix,lsyMatrix-method}
\alias{coerce,lspMatrix,dspMatrix-method}
\alias{coerce,lspMatrix,lgeMatrix-method}
\alias{coerce,lspMatrix,lsyMatrix-method}
\alias{coerce,lsyMatrix,dsyMatrix-method}
\alias{coerce,lsyMatrix,lgeMatrix-method}
\alias{coerce,lsyMatrix,lspMatrix-method}
\alias{coerce,ltCMatrix,lgCMatrix-method}
\alias{coerce,ltCMatrix,ltTMatrix-method}
\alias{coerce,ltTMatrix,dtTMatrix-method}
\alias{coerce,ltTMatrix,lgCMatrix-method}
\alias{coerce,ltTMatrix,lgTMatrix-method}
\alias{coerce,ltTMatrix,ltCMatrix-method}
\alias{coerce,ltTMatrix,ltrMatrix-method}
\alias{coerce,ltpMatrix,dtpMatrix-method}
\alias{coerce,ltpMatrix,lgeMatrix-method}
\alias{coerce,ltpMatrix,ltrMatrix-method}
\alias{coerce,ltrMatrix,dtrMatrix-method}
\alias{coerce,ltrMatrix,lgeMatrix-method}
\alias{coerce,ltrMatrix,ltpMatrix-method}
\alias{coerce,matrix,dgRMatrix-method}
\alias{coerce,matrix,dgTMatrix-method}
\alias{coerce,matrix,dgeMatrix-method}
\alias{coerce,matrix,dsCMatrix-method}
\alias{coerce,matrix,dsTMatrix-method}
\alias{coerce,matrix,dspMatrix-method}
\alias{coerce,matrix,dsyMatrix-method}
\alias{coerce,matrix,dtCMatrix-method}
\alias{coerce,matrix,dtTMatrix-method}
\alias{coerce,matrix,dtpMatrix-method}
\alias{coerce,matrix,dtrMatrix-method}
\alias{coerce,matrix,lgCMatrix-method}
\alias{coerce,matrix,lgTMatrix-method}
\alias{coerce,matrix,lgeMatrix-method}
\alias{coerce,matrix,lsCMatrix-method}
\alias{coerce,matrix,lspMatrix-method}
\alias{coerce,matrix,lsyMatrix-method}
\alias{coerce,matrix,ltCMatrix-method}
\alias{coerce,matrix,ltTMatrix-method}
\alias{coerce,matrix,ltpMatrix-method}
\alias{coerce,matrix,ltrMatrix-method}
\alias{coerce,matrix,ngCMatrix-method}
\alias{coerce,matrix,ngTMatrix-method}
\alias{coerce,matrix,ngeMatrix-method}
\alias{coerce,matrix,nspMatrix-method}
\alias{coerce,matrix,nsyMatrix-method}
\alias{coerce,matrix,ntCMatrix-method}
\alias{coerce,matrix,ntTMatrix-method}
\alias{coerce,matrix,ntpMatrix-method}
\alias{coerce,matrix,ntrMatrix-method}
\alias{coerce,ngCMatrix,dgCMatrix-method}
\alias{coerce,ngCMatrix,lgCMatrix-method}
\alias{coerce,ngCMatrix,ntCMatrix-method}
\alias{coerce,ngTMatrix,dgTMatrix-method}
\alias{coerce,ngTMatrix,lgTMatrix-method}
\alias{coerce,ngTMatrix,lgeMatrix-method}
\alias{coerce,ngTMatrix,ngCMatrix-method}
\alias{coerce,ngTMatrix,ngeMatrix-method}
\alias{coerce,ngTMatrix,ntTMatrix-method}
\alias{coerce,ngeMatrix,dgeMatrix-method}
\alias{coerce,ngeMatrix,lgeMatrix-method}
\alias{coerce,ngeMatrix,ngCMatrix-method}
\alias{coerce,ngeMatrix,ngTMatrix-method}
\alias{coerce,ngeMatrix,nspMatrix-method}
\alias{coerce,ngeMatrix,nsyMatrix-method}
\alias{coerce,ngeMatrix,ntpMatrix-method}
\alias{coerce,ngeMatrix,ntrMatrix-method}
\alias{coerce,nsCMatrix,dsCMatrix-method}
\alias{coerce,nsCMatrix,lsCMatrix-method}
\alias{coerce,nsCMatrix,ngCMatrix-method}
\alias{coerce,nsCMatrix,nsTMatrix-method}
\alias{coerce,nsTMatrix,dsTMatrix-method}
\alias{coerce,nsTMatrix,ngCMatrix-method}
\alias{coerce,nsTMatrix,ngTMatrix-method}
\alias{coerce,nsTMatrix,nsCMatrix-method}
\alias{coerce,nsTMatrix,nsyMatrix-method}
\alias{coerce,nspMatrix,dspMatrix-method}
\alias{coerce,nspMatrix,lspMatrix-method}
\alias{coerce,nspMatrix,ngeMatrix-method}
\alias{coerce,nspMatrix,nsyMatrix-method}
\alias{coerce,nsyMatrix,dsyMatrix-method}
\alias{coerce,nsyMatrix,lsyMatrix-method}
\alias{coerce,nsyMatrix,ngeMatrix-method}
\alias{coerce,nsyMatrix,nspMatrix-method}
\alias{coerce,ntCMatrix,dtCMatrix-method}
\alias{coerce,ntCMatrix,ltCMatrix-method}
\alias{coerce,ntCMatrix,ngCMatrix-method}
\alias{coerce,ntTMatrix,dtTMatrix-method}
\alias{coerce,ntTMatrix,ngCMatrix-method}
\alias{coerce,ntTMatrix,ngTMatrix-method}
\alias{coerce,ntTMatrix,ntCMatrix-method}
\alias{coerce,ntTMatrix,ntrMatrix-method}
\alias{coerce,ntpMatrix,dtpMatrix-method}
\alias{coerce,ntpMatrix,ltpMatrix-method}
\alias{coerce,ntpMatrix,ngeMatrix-method}
\alias{coerce,ntpMatrix,ntrMatrix-method}
\alias{coerce,ntrMatrix,dtrMatrix-method}
\alias{coerce,ntrMatrix,ltrMatrix-method}
\alias{coerce,ntrMatrix,ngeMatrix-method}
\alias{coerce,ntrMatrix,ntpMatrix-method}
\alias{coerce,numLike,dgeMatrix-method}
%
\description{
  These functions are provided for compatibility with older versions
  of \pkg{Matrix} only and may be defunct as soon as the next release.
}
\usage{
%% ## Deprecated in 1.6-0 ... on second thought, perhaps not yet
%% \S4method{chol2inv}{CHMfactor}(x, ...)
%% \S4method{coerce}{CHMfactor,Matrix}(from, to, strict = TRUE)
%% \S4method{coerce}{CHMfactor,CsparseMatrix}(from, to, strict = TRUE)
%% \S4method{coerce}{CHMfactor,RsparseMatrix}(from, to, strict = TRUE)
%% \S4method{coerce}{CHMfactor,TsparseMatrix}(from, to, strict = TRUE)
%% \S4method{coerce}{CHMfactor,dMatrix}(from, to, strict = TRUE)
%% \S4method{coerce}{CHMfactor,dsparseMatrix}(from, to, strict = TRUE)
%% \S4method{coerce}{CHMfactor,pMatrix}(from, to, strict = TRUE)
%% \S4method{coerce}{CHMfactor,sparseMatrix}(from, to, strict = TRUE)
%% \S4method{coerce}{CHMfactor,triangularMatrix}(from, to, strict = TRUE)
%%
## Deprecated in 1.5-4
..2dge(from)
.C2nC(from, isTri)
.T2Cmat(from, isTri)
.asmatrix(x)
.dense2sy(from, ...)
.diag2mat(from)
.diag2sT(from, uplo = "U", kind = ".", drop0 = TRUE)
.diag2tT(from, uplo = "U", kind = ".", drop0 = TRUE)
.dsy2dsp(from)
.dsy2mat(from, keep.dimnames = TRUE)
.dxC2mat(from, chkUdiag)
.m2dgC(from)
.m2lgC(from)
.m2ngC(from)
.m2ngCn(from, na.is.not.0 = FALSE)
.m2ngTn(from, na.is.not.0 = FALSE)
.n2dgT(from)
.nC2d(from)
.nC2l(from)

## Deprecated in 1.5-0
\S4method{coerce}{RsparseMatrix,dgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ddenseMatrix,dgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ddiMatrix,dgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ddiMatrix,dgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ddiMatrix,dtCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgCMatrix,dgTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgCMatrix,dgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgCMatrix,dsCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgCMatrix,dtCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgCMatrix,lgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgCMatrix,ngCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgTMatrix,dgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgTMatrix,dgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgTMatrix,dsTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgTMatrix,dtCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgTMatrix,dtTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgeMatrix,dgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgeMatrix,dgTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgeMatrix,dsTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgeMatrix,dspMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgeMatrix,dsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgeMatrix,dtrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dgeMatrix,lgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsCMatrix,dgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsCMatrix,dgTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsCMatrix,dgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsCMatrix,dsRMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsCMatrix,dsTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsCMatrix,dsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsCMatrix,lsCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsCMatrix,nsCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsTMatrix,dgTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsTMatrix,dgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsTMatrix,dsCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsTMatrix,dsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsTMatrix,lsTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dspMatrix,dsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dspMatrix,lspMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsyMatrix,dsCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsyMatrix,dsTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsyMatrix,dspMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dsyMatrix,lsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtCMatrix,dgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtCMatrix,dgTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtCMatrix,dgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtCMatrix,dsCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtCMatrix,dtTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtCMatrix,dtrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtCMatrix,ltCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtCMatrix,ntCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtTMatrix,dgTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtTMatrix,dgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtTMatrix,dtCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtTMatrix,dtrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtpMatrix,dtTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtpMatrix,dtrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtpMatrix,ltpMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtrMatrix,dtpMatrix}(from, to, strict = TRUE)
\S4method{coerce}{dtrMatrix,ltrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{indMatrix,ngTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{indMatrix,ngeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lMatrix,dgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lgCMatrix,dgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lgCMatrix,lgTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lgCMatrix,lgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lgCMatrix,ltCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lgTMatrix,dgTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lgTMatrix,lgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lgTMatrix,lgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lgTMatrix,lsCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lgTMatrix,ltTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lgeMatrix,dgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lgeMatrix,lgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lgeMatrix,lgTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lgeMatrix,lspMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lgeMatrix,lsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lgeMatrix,ltpMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lgeMatrix,ltrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lsCMatrix,dsCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lsCMatrix,lgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lsCMatrix,lgTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lsCMatrix,lsTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lsTMatrix,lgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lsTMatrix,lgTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lsTMatrix,lsCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lsTMatrix,lsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lspMatrix,dspMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lspMatrix,lgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lspMatrix,lsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lsyMatrix,dsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lsyMatrix,lgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{lsyMatrix,lspMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ltCMatrix,lgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ltCMatrix,ltTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ltTMatrix,dtTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ltTMatrix,lgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ltTMatrix,lgTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ltTMatrix,ltCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ltTMatrix,ltrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ltpMatrix,dtpMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ltpMatrix,lgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ltpMatrix,ltrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ltrMatrix,dtrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ltrMatrix,lgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ltrMatrix,ltpMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,dgRMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,dgTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,dgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,dsCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,dsTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,dspMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,dsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,dtCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,dtTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,dtpMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,dtrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,lgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,lgTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,lgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,lsCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,lspMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,lsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,ltCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,ltTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,ltpMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,ltrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,ngCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,ngTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,ngeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,nspMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,nsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,ntCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,ntTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,ntpMatrix}(from, to, strict = TRUE)
\S4method{coerce}{matrix,ntrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngCMatrix,dgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngCMatrix,lgCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngCMatrix,ntCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngTMatrix,dgTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngTMatrix,lgTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngTMatrix,lgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngTMatrix,ngCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngTMatrix,ngeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngTMatrix,ntTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngeMatrix,dgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngeMatrix,lgeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngeMatrix,ngCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngeMatrix,ngTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngeMatrix,nspMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngeMatrix,nsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngeMatrix,ntpMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ngeMatrix,ntrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nsCMatrix,dsCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nsCMatrix,lsCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nsCMatrix,ngCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nsCMatrix,nsTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nsTMatrix,dsTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nsTMatrix,ngCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nsTMatrix,ngTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nsTMatrix,nsCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nsTMatrix,nsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nspMatrix,dspMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nspMatrix,lspMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nspMatrix,ngeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nspMatrix,nsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nsyMatrix,dsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nsyMatrix,lsyMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nsyMatrix,ngeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{nsyMatrix,nspMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ntCMatrix,dtCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ntCMatrix,ltCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ntCMatrix,ngCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ntTMatrix,dtTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ntTMatrix,ngCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ntTMatrix,ngTMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ntTMatrix,ntCMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ntTMatrix,ntrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ntpMatrix,dtpMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ntpMatrix,ltpMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ntpMatrix,ngeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ntpMatrix,ntrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ntrMatrix,dtrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ntrMatrix,ltrMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ntrMatrix,ngeMatrix}(from, to, strict = TRUE)
\S4method{coerce}{ntrMatrix,ntpMatrix}(from, to, strict = TRUE)
\S4method{coerce}{numLike,dgeMatrix}(from, to, strict = TRUE)
}
\seealso{
  \code{\link{Deprecated}},
  \code{\link{base-deprecated}},
  \code{\link{Matrix-defunct}}
}
