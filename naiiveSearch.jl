import StaticArrays:MMatrix
import Combinatorics:multiset_permutations
import Printf:@printf
function main()
    A=zeros(MMatrix{4,4,Int8})
    max=0.0
    maxMat=deepcopy(A)
    for a in multiset_permutations([10,9,8,7,6,5,4,3,2,1],[5,3,3,3,2,2,3,3,3,4],10)
        A[1,1]=a[1]
        A[2,1]=a[2]
        A[3,1]=a[3]
        A[4,1]=a[4]
        A[1,2]=a[2]
        A[2,2]=a[5]
        A[3,2]=a[6]
        A[4,2]=a[7]
        A[1,3]=a[3]
        A[2,3]=a[6]
        A[3,3]=a[8]
        A[4,3]=a[9]
        A[1,4]=a[4]
        A[2,4]=a[7]
        A[3,4]=a[9]
        A[4,4]=a[10]
        if abs(det(A))<0.5
            continue
        end
        condA=cond(A)
        if condA>max
            max=condA
            maxMat=deepcopy(A)
            @printf("%s, %f\n",A,max)
        end
    end
end
