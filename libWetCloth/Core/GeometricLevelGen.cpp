//
// This file is part of the libWetCloth open source project
//
// The code is licensed under the same terms as a Clear BSD License but further
// restricted to academic and non-commercial use (commercial licenses may be
// obtained by contacting the faculty of the Columbia Computer Graphics Group
// or Columbia Technology Ventures).
//
// Copyright 2015 Xinxin Zhang
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
// to endorse or promote products derived from this software without specific
// prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
// LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

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

	int nni = ceil((float)ni/2.0);
	int nnj = ceil((float)nj/2.0);
	int nnk = ceil((float)nk/2.0);
	SparseMatrix<T> r;
	SparseMatrix<T> p;
	p.resize(ni*nj*nk);
	p.zero();
	r.resize(nni*nnj*nnk);
	r.zero();

	for(int k=0;k<nnk;k++)for(int j=0;j<nnj;j++)for (int i=0;i<nni;i++)
	{
		unsigned int index = (k*nnj +j)*nni + i;
		for(int kk=0;kk<=1;kk++)for(int jj=0;jj<=1;jj++)for(int ii=0;ii<=1;ii++)
		{
			int iii = i*2 + ii;
			int jjj = j*2 + jj;
			int kkk = k*2 + kk;
			if(iii<ni && jjj<nj && kkk<nk)
			{
				unsigned int index2 = (kkk*nj + jjj)*ni + iii;
				r.set_element(index, index2, (T)0.125);
				p.set_element(index2,index, 1.0);
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
