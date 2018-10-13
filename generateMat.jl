import LinearAlgebra:cond,reshape,Matrix,UpperTriangular,diag,I,det
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

#function naiiveMatIter()
function randMat4(n)
		return MMatrix{4,4}([mod(i,n)+1 for i in rand(Int,(4,4))])
end

function randMat4_NS(n)
		A=randMat4(n)
		while abs(mod(det(A),n))<0.1
				A=randMat4(n)
		end
		return A
end

function randSMat4(n)
		A=randMat4(n)
		return MMatrix{4,4}([mod(i,n)+1 for i in A+transpose(A)-diag(A).*Matrix(1*I,4,4)])
end

function randSMat4_NS(n)
		A=randSMat4(n)
		while abs(mod(det(A),n))<0.1
				A=randSMat4(n)
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
detB^(i+n*j)*detA === 1 (mod p), n in N+
"""
function solveModEqn(detA::Int,detB::Int,p::Int)::Int

		i=1
		j=1
		while mod(detB^(2*i)*detA,p)!=1
				i+=1
		end
		while mod(detB^(2*j),p)!=1
				j+=1
		end
		return i,j
end

function matModExp(A::MMatrix,e::Int,n::Int)
		Aexp=deepcopy(A)
		for i=1:(e-1)
				Aexp=mod.(Aexp*A,n)
		end
		return Aexp
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
	ii,jj=solveModEqn(Int(det(A)),Int(det(B)),n)
	A=ModCong!(A,matModExp(B,ii,n),n)
	B=matModExp(B,jj,n)
	Aorigin=deepcopy(A)
	#B=deepcopy(A)
	max=0.0
	maxMat=zeros(4)
	for i=1:N
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
		max=0.0
		maxi=0.0
		maxMat=zeros(Int,4)
		maxMati=zeros(Int,4)
		B=zeros(Int,4)
		Bi=zeros(Int,4)
		A=randSMat4_NS(10)
		for i=1:N
				if mod(i,100)==0
						A=randSMat4_NS(10)
				end
				if i-floor(sqrt(i))^2<1
						@printf("Current i=%d,max=%f\n",i,max)
				end
				#h=rand01TriMat()
				maxi,maxMati,Bi=gensubgroup(100000,A,randMat4_NS(10))
				if maxi>max
						max=maxi
						maxMat=maxMati
						B=Bi
				end
		end
		println(max)
		println(maxMat)
		println(B)
		@printf("Det(maxMat)=%f,Det(B)=%f,cond(B)=%f\n",det(maxMat),det(B),cond(B))
end
