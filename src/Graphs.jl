include("Test.jl")
using Plots, LinearAlgebra, LaTeXStrings

#Gráfico de dispersão para erro relativo com um dado método
function scatter_error(x,y,metodo)
    scatter(x,y,label=metodo)
    xlabel!(L"$n$")
    ylabel!(L"$\frac{||x-x^*||}{||x^*||}$")
    savefig("images/Erro Relativo - "*metodo*".png")
end

#Gráficos de dispersão para resíduo relativo
function scatter_residue(x,y,metodo)
    scatter(x,y,label=metodo)
    xlabel!(L"$n$")
    ylabel!(L"$\frac{||Ax-b||}{||b||}$")
    savefig("images/Residuos - "*metodo*".png")
end

#Gráficos de dispersão para tempo
function scatter_time(x,y,metodo)
    scatter(x,y,label=metodo)
    xlabel!(L"$n$")
    ylabel!(L"$Tempo (s)$")
    savefig("images/Tempos - "*metodo*".png")
end

#Gráficos de reta/dispersão para um dado método com o desempenho de tempo do algoritmo, residuo e erro relativo para um conjunto de dimensões n
function plot_scatter(dim)
    #Eixo x do gráfico são as dimensões
    n=size(dim,1)

    #Pegando os resultados
    resultPLU, resultJac=runTest(dim)

    #Plotando e salvando os gráficos dos erros
    scatter_error(dim,resultJac[:,1],"Gauss-Jacobi")
    scatter_residue(dim,resultJac[:,2],"Gauss-Jacobi")
    scatter_time(dim,resultJac[:,3],"Gauss-Jacobi")
    scatter_error(dim,resultPLU[:,1],"LU com Piv. Parcial")
    scatter_residue(dim,resultPLU[:,2],"LU com Piv. Parcial")
    scatter_time(dim,resultPLU[:,3],"LU com Piv. Parcial")
end

dim=[i for i=1000:1000:5000]
#plot_scatter(dim)
