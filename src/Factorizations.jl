using LinearAlgebra
#Recebe a matriz A, retorna L,U e matriz de permutação P -> Implementação da LU da questão 1)
function LUP(A)
    m=size(A,1)
    L=Matrix{Float64}(I,m,m)
    P=Matrix{Float64}(I,m,m)
    U=copy(A);
    for j=1:m-1
        #Determinando o pivo
        ind=argmax(broadcast(abs,U[j:m,j]))+j-1
        #Permutando
        U[[j,ind],:]=U[[ind,j],:]
        L[[j,ind],1:j-1]=L[[ind,j],1:j-1]
        P[[j,ind],:]=P[[ind,j],:]
        #Eliminando
        for i=j+1:m
            L[i,j]=U[i,j]/U[j,j]
            U[i,:]=U[i,:]-L[i,j]*U[j,:]
        end
    end
    return L,U,P
end

#Recebe a matriz A, retorna L,U e matriz de permutação P
#Versão do ChatGPT -> Implementação para a questão 2
function LUPGPT(A)
    n, m = size(A)
    if n != m
        error("Input matrix must be square")
    end
    L = Matrix{Float64}(I, n, n)  # Initialize L matrix as identity matrix
    U = copy(A)                    # Initialize U matrix as a copy of the input matrix A
    P = Matrix{Float64}(I, n, n)  # Initialize P matrix as identity matrix
    for k = 1:n-1
        pivot_row = argmax(abs.(U[k:n, k]))[1] + k - 1  # Find row with maximum absolute value in the column
        if pivot_row != k
            # Swap rows in U matrix
            U[[k, pivot_row], k:end] = U[[pivot_row, k], k:end]

            # Swap rows in P matrix
            P[[k, pivot_row], :] = P[[pivot_row, k], :]
            
            if k > 1
                # Swap rows in L matrix only if k > 1
                L[[k, pivot_row], 1:k-1] = L[[pivot_row, k], 1:k-1]
            end
        end

        for i = k+1:n
            factor = U[i, k] / U[k, k]
            L[i, k] = factor  # Store the factor in L matrix
            U[i, k:end] -= factor * U[k, k:end]
        end
    end
    return L, U, P
end
