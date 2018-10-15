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
function ModConj!(A::MMatrix{4,4,Int},B::MMatrix{4,4,Int},invB::SMatrix{4,4,Int},n::Int)::MMatrix
	temp=invB*A*B
	for i=1:16
		@inbounds A[i]=mod(temp[i],n)
	end
	#return [mod(i,n) for i in transpose(B)*A*B]
	return A
end

function ModMult!(A::MMatrix{4,4,Int},B::MMatrix{4,4,Int},n::Int)::MMatrix
    temp=A*B
    for i=1:16
        @inbounds A[i]=mod(temp[i],n)
    end
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
    moddetA=mod(detA,p)
    moddetB2=mod(mod(detB,p)^2,p)
		pq=mod(moddetB2*moddetA,p)
		q2=moddetB2
		while (q2!=1)&(q2!=mod(-1,p))
				q2=mod(q2*moddetB2,p)
				i+=1
		end
		while (pq!=1)&(pq!=mod(-1,p))
				pq=mod(pq*moddetB2,p)
				j+=1
				if j>i
					j=-1
					break
				end
		end
		return i,j
end
"""
detA^i == +-1 (mod p)
"""
function solveModEqn(detA::Int,p::Int)
    i=1
    moddetA=mod(detA,p)
    detAn=moddetA
    while (detAn!=1)&(detAn!=-1)
        detAn=mod(detAn*moddetA,p)
        i+=1
    end
    return i
end


function matModExp(A::MMatrix,e::Int,n::Int)::MMatrix
		Aexp=deepcopy(A)
		for i=1:(e-1)
				ModMult!(Aexp,A,n)
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
  #@assert (Int(mod(det(A),11))==1) & (Int(mod(det(A),11))==1)
	n=11
	Aorigin=deepcopy(A)
	#invB=SMatrix{4,4,Int}(round.(inv(B)))
	max=0.0
	maxMat=zeros(4)
	for i=1:N
		#ModConj!(A,B,invB,n)
		ModCong!(A,B,n)
		#if abs(det(A))<0.9
			#print("Singular Matrix i=")
			#println(i)
		#	break
		#end
    condA=cond(A)
		if condA>max
			max=condA
			#println(A)
			maxMat=deepcopy(A)
		end
		if A==Aorigin
			#print("Cycled i=")
			#println(i)
			break
		end
	end
  return max,maxMat,B
end
"""
det(A)===1 (mod 11)
A->A*A
"""
function gensubgroup(N,A::MMatrix)
    #@assert Int(mod(det(A),11))==1
    n=11
    Aorigin=deepcopy(A)
    max=0.0
    maxMat=zeros(MMatrix{4,4,Int})
    for i=1:N
        ModMult!(A,Aorigin,n)
        condA=cond(A)
        if condA>max
            max=condA
            maxMat=deepcopy(A)
        end
        if A==Aorigin
            #@printf("Cycled i=%d\n",i)
            break
        end
    end
    return max,maxMat
end


"""
Finding 4x4 Matrix with elements of 0-10 that have big condition number
"""
function genBigCond(N)
    n=11
    max=0.0
    maxi=0.0
    maxMat=zeros(MMatrix{4,4,Int})
    maxMati=zeros(MMatrix{4,4,Int})
    maxB=zeros(MMatrix{4,4,Int})
    Bi=zeros(MMatrix{4,4,Int})
		f=open("logBigCond_NonS.txt","a")
    for i=1:N
        A=randMat4_NS(n)
        ii=solveModEqn(Int(det(A)),11)
        A=matModExp(A,ii,n)
        maxi,maxMati=gensubgroup(100000,A)
        if maxi>max
            max=maxi
            maxMat=maxMati
            if max>0
                @printf("i=%d/%d,\n%s %f,det(A)=%f\n",i,N,maxMat,max,det(maxMat))
                @printf(f,"%d %s %f\n",i,maxMat,max)
                flush(f)
            end
        end
    end
    close(f)
    @printf("det(maxMat)=%f",det(maxMat))
end

"""
Finding *Symmetric* 4x4 Matrix with elements of 0-10 that have big condition number
"""
function genBigCondS(N)
		n=11
		max=0.0
		maxi=0.0
		maxMat=zeros(Int,4)
		maxMati=zeros(Int,4)
		maxB=zeros(Int,4)
		Bi=zeros(Int,4)
		A=randSMat4_NS(n)
		f=open("logBigCond.txt","a")
		for i=1:N
				if mod(i,100)==0
						A=randSMat4_NS(n)
				end
				#if i-floor(i^(1/2))^2<1
				#		@printf("Current i=%d/%d,max=%f\n",i,N,max)
				#end
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
							@printf("i=%d/%d,\n%s %f\n %s\n",i,N,maxMat,max,maxB)
							@printf(f,"%d %s %f %s\n",i,maxMat,max,maxB)
							flush(f)
						end
				end
		end
		close(f)
		println(max)
		println(maxMat)
		println(maxB)
		@printf("Det(maxMat)=%f,Det(B)=%f,cond(B)=%f\n",det(maxMat),det(maxB),cond(maxB))
end
