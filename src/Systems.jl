using LinearAlgebra
include("Factorizations.jl")

#Gera a matriz de teste da questão 3)
function testMatrix(n)
    A=zeros(n,n)
    for i=1:n
        A[i,i]=10
        if i>1
            A[i,1]=1/i
            A[1,i]=1/n
        end
    end
    return A
end

#Vetor de teste b para o sistema linear - questão 4)
function testb(n)
    A=testMatrix(n); x=2*rand(n,1).-1; b=A*x
    return b,x
end

#Matriz só tem entradas na primeira linha, primeira coluna e na diagonal - Método adaptado da 4.A)
#Podemos só referenciar esses elementos em cada iteração, em vez de ter que trabalhar com a inversa de D etc.
#Função resolve o sistema da questão pelo método de jacobi adaptado
function jacobiAdaptado(x0,b,maxIt=1000,eps=1e-8)
    n=size(x0,1); x=copy(x0); xNew=copy(x); dx=0; it=0;
    while it<maxIt
        xNew[1]=(b[1]-(sum(x[2:n]))/n)/10
        for j=2:n
            xNew[j]=(b[j]-x[1]/j)/10
        end
        dx=norm(x-xNew)
        x.=xNew; it+=1
        if dx<eps
            break
        end
    end
    A=testMatrix(n)
    return x
end

#Método adaptado feito pelo chatgpt - 4b) -> faz a mesma coisa que a função anterior.
function jacobiAdaptadoGPT(x0, b, max_iter=1000, tol=1e-8)
    n = length(x0)
    x = copy(x0)
    x_new = zeros(n)
    iter_count = 0
    for k in 1:max_iter
        # Update x[1]
        sum1 = sum(x[2:end])  # Sum of elements from x[2] to x[n]
        x_new[1] = (b[1] - sum1 / n) / 10
        # Update x[i] for i > 1
        for j in 2:n
            x_new[j] = (b[j] / 10 - 1 / (10 * j) * x[1])
        end
        dx = norm(x_new - x)
        if dx < tol
            return x_new
        end
        x .= x_new
        iter_count += 1
    end
    return x
end

#Função que resolve um sistema Ax=B usando LU com pivoteamento parcial
#opt=1 -> Usa a fatoração LU feita, caso contrário, usa a do chatgpt
function resolveLUP(A,b,opt=1)
    n=size(b,1)
    #Vendo qual fatoração usar
    if opt==1
        L,U,P=LUP(A)
    else
        L,U,P=LUPGPT(A)
    end
    #Variáveis e vetores prealocados
    b_aux=P*b; soma=0; x=zeros(n,1); y=zeros(n,1)
    #Resolvemos os dois sistemas triangulares
    #Inferior Ly=b_aux, y=Ux
    for i=1:n
        for j=1:(i-1)
            soma+=L[i,j]/L[i,i]*y[j]
        end
        y[i]=b_aux[i]/L[i,i]-soma
        soma=0
    end
    #Superior Ux=y
    for i=n:-1:1
        for j=i+1:n
            soma+=U[i,j]/U[i,i]*x[j]
        end
        x[i]=y[i]/U[i,i]-soma
        soma=0
    end
    return x
end