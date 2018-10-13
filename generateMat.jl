import LinearAlgebra:cond,reshape,Matrix,UpperTriangular,diag,I,det,inv
import Printf:@printf
import StaticArrays:SMatrix,MMatrix
function ModCong!(A::MMatrix{4,4,Int},B::MMatrix{4,4,Int},n::Int)::MMatrix
	temp=transpose(B)*A*B
	for i=1:16
		@inbounds A[i]=mod(temp[i],n)
	end
	#return [mod(i,n) for i in transpose(B)*A*B]
	return A
end
"""
Convergent rate of Conj is slower than Cong
"""
function ModConj!(A::MMatrix,B::MMatrix,invB::SMatrix{4,4,Int},n::Int)::MMatrix
	temp=invB*A*B
	for i=1:16
		@inbounds A[i]=mod(temp[i],n)
	end
	#return [mod(i,n) for i in transpose(B)*A*B]
	return A
end

#function naiiveMatIter()
function randMat4(n)
		return MMatrix{4,4}([mod(i,n)+1 for i in rand(Int,(4,4))])
end

function randMat4_NS(n)
		A=randMat4(n-1)
		while abs(mod(det(A),n))<0.1
				A=randMat4(n-1)
		end
		return A
end

function randSMat4(n)
		A=randMat4(n)
		return MMatrix{4,4}([mod(i,n)+1 for i in A+transpose(A)-diag(A).*Matrix(1*I,4,4)])
end

function randSMat4_NS(n)
		A=randSMat4(n-1)
		while abs(mod(det(A),n))<0.1
				A=randSMat4(n-1)
		end
		return A
end

function randSMat4_Det1(n)
		A=randSMat4_NS(n)
		Aorigin=deepcopy(A)
		while abs(det(A)-1)>0.5
				ModCong!(A,Aorigin,n)
		end
		return A
end

function rand01TriMat()::MMatrix
		A=[mod(i,2) for i in rand(Int,(4,4))]
		return MMatrix{4,4}(UpperTriangular(A-diag(A).*Matrix(1*I,4,4)+Matrix(1*I,4,4)))
end

"""
detB^(j+n*i)*detA === 1 (mod p), n in N+
"""
function solveModEqn(detA::Int,detB::Int,p::Int)
		i=1
		j=1
		pq=mod(mod(detB,p)^2*mod(detA,p),p)
		q2=mod(mod(detB,p)^2,p)
		while q2!=1
				q2=mod(q2*mod(detB,p)^2,p)
				i+=1
		end
		while pq!=1
				pq=mod(pq*mod(detB,p)^2,p)
				j+=1
				if j>i
					j=-1
					break
				end
		end
		return i,j
end

function matModExp(A::MMatrix,e::Int,n::Int)::MMatrix
		Aexp=deepcopy(A)
		for i=1:(e-1)
				Aexp=mod.(Aexp*A,n)
		end
		return MMatrix{4,4}(Aexp)
end

"""
det(A)===1 (mod 11)
det(B)===1 (mod 11)
A->B^T*A*B
"""
function gensubgroup(N,A::MMatrix,B::MMatrix)
	#h=rand01TriMat()
	#h=[1 0 1 0;0 1 0 0;0 0 1 0;0 0 0 1]
	#B=[1 3 10 10;3 4 8 9;10 8 3 9;10 9 9 3]
	#A=MMatrix{4,4}(transpose(h)*h)#[2 7 10 10;7 10 10 9;10 10 10 1;10 9 1 9]
	#A=randSMat4_NS(11)
	n=11
	Aorigin=deepcopy(A)
	#invB=SMatrix{4,4,Int}(round.(inv(B)))
	#B=deepcopy(A)
	max=0.0
	maxMat=zeros(4)
	for i=1:N
		#ModConj!(A,B,invB,n)
		ModCong!(A,B,n)
		if abs(det(A))<0.9
			#print("Singular Matrix i=")
			#println(i)
			break
		end
		if cond(A)>max
			max=cond(A)
			#println(A)
			maxMat=deepcopy(A)
		end
		if A==Aorigin
			#print("Cycled i=")
			#println(i)
			break
		end
	end
	#println(max)
	#println(maxMat)
	#println(B)
  return max,maxMat,B
end

function main(N)
		n=11
		max=0.0
		maxi=0.0
		maxMat=zeros(Int,4)
		maxMati=zeros(Int,4)
		maxB=zeros(Int,4)
		Bi=zeros(Int,4)
		A=randSMat4_NS(n)
		for i=1:N
				if mod(i,100)==0
						A=randSMat4_NS(n)
				end
				if i-floor(i^(1/2))^2<1
						@printf("Current i=%d/%d,max=%f\n",i,N,max)
				end
				#h=rand01TriMat()
				B=randMat4_NS(n)
				ii,jj=solveModEqn(Int(det(A)),Int(det(B)),n)
				while ((mod(det(B),n)==1)&(mod(det(A),n)!=1))|(jj==-1)
					A=randSMat4_NS(n)
					ii,jj=solveModEqn(Int(det(A)),Int(det(B)),n)
				end
				#println(mod(det(A),n),", ",mod(det(B),n))
				A=ModCong!(A,matModExp(B,jj,n),n)
				B=matModExp(B,2*ii,n)
				#println(mod(det(A),n),", ",mod(det(B),n),"////")
				maxi,maxMati,Bi=gensubgroup(100000,A,B)
				if maxi>max
						max=maxi
						maxMat=maxMati
						maxB=Bi
						if max>70000
							println(maxMat)
						end
				end
		end
		println(max)
		println(maxMat)
		println(maxB)
		@printf("Det(maxMat)=%f,Det(B)=%f,cond(B)=%f\n",det(maxMat),det(maxB),cond(maxB))
end
