//
// This file is part of the libWetCloth open source project
//
// The code is licensed solely for academic and non-commercial use under the
// terms of the Clear BSD License. The terms of the Clear BSD License are
// provided below. Other licenses may be obtained by contacting the faculty
// of the Columbia Computer Graphics Group or a Columbia University licensing officer.
//
// We would like to hear from you if you appreciate this work.
//
// The Clear BSD License
//
// Copyright 2015 Xinxin Zhang
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
//  to endorse or promote products derived from this software without specific
//  prior written permission.
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


template<class T>
void levelGen<T>::generateLevelsGalerkinCoarsening(vector<FixedSparseMatrix<T> > &A_L,
	vector<FixedSparseMatrix<T> > &R_L,
	vector<FixedSparseMatrix<T> > &P_L,
	vector<Vector3i> &S_L, int & total_level, 
	FixedSparseMatrix<T> &A, 
	int ni,int nj,int nk)
{
#ifdef AMG_VERBOSE
	cout<<"building levels ...... "<<endl;
#endif
	A_L.resize(0);
	R_L.resize(0);
	P_L.resize(0);
	S_L.resize(0);
	total_level = 1;
  A_L.push_back(std::move(A));
	S_L.push_back(Vector3i(ni,nj,nk));
	int nni = ni, nnj = nj, nnk = nk;
	unsigned int unknowns = ni*nj*nk;
	while (unknowns > 16*16*16)
	{
		A_L.push_back(FixedSparseMatrix<T>());
		R_L.push_back(FixedSparseMatrix<T>());
		P_L.push_back(FixedSparseMatrix<T>());
		nni = ceil((float)nni/2.0);
		nnj = ceil((float)nnj/2.0);
		nnk = ceil((float)nnk/2.0);

		S_L.push_back(Vector3i(nni,nnj,nnk));
		unknowns = nni*nnj*nnk;
		total_level++;
	}
	
	for (int i=0;i<total_level-1;i++)
	{
		generateRP((A_L[i]), (R_L[i]),(P_L[i]),S_L[i][0], S_L[i][1],S_L[i][2]);
		FixedSparseMatrix<T> temp;
		multiplyMat((A_L[i]),(P_L[i]),temp,1.0);
		multiplyMat((R_L[i]),temp, (A_L[i+1]),0.5);
		temp.resize(0);
		temp.clear();
	}
#ifdef AMG_VERBOSE
	cout<<"build levels done"<<endl;
#endif
}

template struct levelGen<double>;
