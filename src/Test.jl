include("Factorizations.jl")
include("Systems.jl")

#Função para rodar os testes da questão 5)
function runTest(n_conj)
    n=size(n_conj,1)
    #Matriz com os erros, resíduos e tempo de execução para cada método
    resultPLU=zeros(n,3); resultJac=zeros(n,3)

    #Rodando os testes
    for i=1:n
        #Criando o sistema Ax=b
        A=testMatrix(n_conj[i])
        b,x_opt=testb(n_conj[i])

        #Vetor auxiliar para armazenar as respostas obtidas
        x=zeros(n_conj[i],1)
        t=0; err=0; res=0
        #Origem é o ponto inicial
        x0=zeros(n_conj[i],1)

        #Resolvendo com LUP
        t=@elapsed x=resolveLUP(A,b,1) #LUP
        err=norm(x-x_opt)/norm(x_opt)
        res=norm(A*x-b)/norm(b)
        
        #Colocando na matriz
        resultPLU[i,:]=[err,res,t]

        #resolvendo com Jacobi
        t=@elapsed x=jacobiAdaptado(x0,b)
        err=norm(x-x_opt)/norm(x_opt) 
        res=norm(A*x-b)/norm(b)

        #Colocando na matriz
        resultJac[i,:]=[err,res,t] 
    end
    return resultPLU, resultJac
end