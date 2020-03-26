//
// This file is part of the libWetCloth open source project
//
// Copyright 2015 Xinxin Zhang
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "GeometricLevelGen.h"


//#define AMG_VERBOSE

template<class T>
void levelGen<T>::generateRP(const FixedSparseMatrix<T> &A,
                             FixedSparseMatrix<T> &R,
                             FixedSparseMatrix<T> &P,
                             int ni, int nj, int nk)
{
	//matlab code(2D):
	//m = ceil(M/2);
	//n = ceil(N/2);
	//R = sparse(zeros(m*n, M*N));
	//coef = [1 3 3 1; 3 9 9 3; 3 9 9 3; 1 3 3 1]/16.0;
	//for i=0:m-1
	//	for j=0:n-1
	//		index = i*n + j +1;
	//		for ii=0:1
	//			for jj=0:1
	//				iii = i*2+ii;
	//				jjj = j*2+jj;
	//				if iii<M && jjj<N && iii>=0 && jjj>=0
	//					index2 = iii*N + jjj +1;
	//					if A(index2,index2)~=0
	//						R(index, index2) = 0.25;
	//						%R(index, index2) = coef(ii+2,jj+2)/4;
	//					else
	//						R(index, index2) = 0.0;
	//					end
	//				end
	//			end
	//		end
	//	end
	//end

	//P = 4*R';

	int nni = ceil((float)ni / 2.0);
	int nnj = ceil((float)nj / 2.0);
	int nnk = ceil((float)nk / 2.0);
	SparseMatrix<T> r;
	SparseMatrix<T> p;
	p.resize(ni * nj * nk);
	p.zero();
	r.resize(nni * nnj * nnk);
	r.zero();

	for (int k = 0; k < nnk; k++)for (int j = 0; j < nnj; j++)for (int i = 0; i < nni; i++)
			{
				unsigned int index = (k * nnj + j) * nni + i;
				for (int kk = 0; kk <= 1; kk++)for (int jj = 0; jj <= 1; jj++)for (int ii = 0; ii <= 1; ii++)
						{
							int iii = i * 2 + ii;
							int jjj = j * 2 + jj;
							int kkk = k * 2 + kk;
							if (iii < ni && jjj < nj && kkk < nk)
							{
								unsigned int index2 = (kkk * nj + jjj) * ni + iii;
								r.set_element(index, index2, (T)0.125);
								p.set_element(index2, index, 1.0);
							}
						}
			}

	R.construct_from_matrix(r);
	P.construct_from_matrix(p);
	r.clear();
	p.clear();

	//transposeMat(R,P,(T)8.0);


}


template struct levelGen<double>;
