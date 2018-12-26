module NURBS

using IntervalSets
using Luxor
using ElementaryCalculus
import ParametricDraw.ChangeUnit
import ParametricDraw.BézPts
import ParametricDraw.LxrPt

export NURBS1mfd, NURBS2mfd, Bs, Ḃs, Bsupp, NURBSmapping, href, pref, NURBSdraw, N, N′, NURBSṁapping

mutable struct NURBS1mfd
    p::Int64
    k::Array{Float64,1}
    w::Array{Float64,1}
    a::Array{Float64,2}
    function NURBS1mfd(p,k,w,a)
        if (size(a)≠(length(k)-p-1,2))
            error("dim-error")
        elseif (size(w)≠(length(k)-p-1,))
            error("dim-error2")
        elseif (!issorted(k))
            error("knots not sorted")
        else
            new(p,k,w,a)
        end
    end
end

mutable struct NURBS2mfd
    p::Array{Int64,1}
    k::Array{Array{Float64,1},1}
    w::Array{Float64,2}
    a::Array{Float64,3}
    function NURBS2mfd(p,k,w,a)
        if (length(p)≠2)
            error("p-error")
        elseif (length(k)≠2)
            error("k-error")
        elseif ([size(w)...]≠(length.(k)-p.-1))
            error("dim-error")
        elseif (size(a)≠((length.(k)-p.-1)...,2))
            error("dim-error")
        elseif (!*(issorted.(k)...))
            error("knots not sorted")
        else
            new(p,k,w,a)
        end
    end
end

function Bs(i::Int64,p::Int64,k,t)::Float64
    if(p==0)
        return k[i]≤t<k[i+1]||(k[i]≠k[i+1]==k[end]==t)
    else
        return (((k[i+p]-k[i]≠0) ? Bs(i,p-1,k,t)*(t-k[i])/(k[i+p]-k[i]) : 0)
        +((k[i+p+1]-k[i+1]≠0) ? Bs(i+1,p-1,k,t)*(k[i+p+1]-t)/(k[i+p+1]-k[i+1]) : 0))
    end
end

function Ḃs(i::Int64,p::Int64,k,t)::Float64
    return p*(((k[i+p]-k[i]≠0) ? Bs(i,p-1,k,t)/(k[i+p]-k[i]) : 0)
    -((k[i+p+1]-k[i+1]≠0) ? Bs(i+1,p-1,k,t)/(k[i+p+1]-k[i+1]) : 0))
end

function Bsupp(i,p,k)::ClosedInterval
    return k[i]..k[i+p+1]
end

function BsCoef(f,p::Int64,k::Array{T,1}) where T<:Real
    n=length(k)-p-1
    κ=[((n+1-i)*k[i]+i*k[i+p+1])/(n+1) for i ∈ 1:n]
    [Bs(j,p,k,κ[i]) for i ∈ 1:n, j ∈ 1:n]\f.(κ)
end

function NURBSmapping(N1::NURBS1mfd,t)
    p=N1.p
    k=N1.k
    w=N1.w
    a=N1.a
    n=length(k)-p-1
    return sum(Bs(I,p,k,t)*a[I,:]*w[I] for I ∈ 1:n)/sum(Bs(J,p,k,t)*w[J] for J ∈ 1:n)
end

function NURBSmapping(N2::NURBS2mfd,u)
    p₁,p₂=p=N2.p
    k₁,k₂=k=N2.k
    w=N2.w
    a=N2.a
    n₁,n₂=n=length.(k)-p.-1
    return sum(Bs(I₁,p₁,k₁,u[1])*Bs(I₂,p₂,k₂,u[2])*a[I₁,I₂,:]*w[I₁,I₂] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂)/sum(Bs(J₁,p₁,k₁,u[1])*Bs(J₂,p₂,k₂,u[2])*w[J₁,J₂] for J₁ ∈ 1:n₁, J₂ ∈ 1:n₂)
end

function NURBSdraw(N1::NURBS1mfd;filename="NURBS1mfd.svg",up=5,down=-5,right=5,left=-5,zoom=1,unitlength=(100,"pt"),points=true,partitionnumber=5)
    step, unit=(unitlength[1]*zoom,unitlength[2])
    Drawing(step*(right-left),step*(up-down),filename)

    Luxor.origin(-step*left,step*up)
    setline(zoom)
    background("white")

    p=N1.p
    k=N1.k
    a=N1.a
    n=length(k)-p-1
    𝒑(t)=NURBSmapping(N1,t)


    K=DelDpl(k[1+p:end-p])
    KK=Float64[]
    for i in 1:length(K)-1
        for j in 0:partitionnumber
            push!(KK,((partitionnumber-j)*K[i]+j*K[i+1])/partitionnumber)
        end
    end
    K=DelDpl(KK)
    N=length(K)-1

    sethue("red")
    drawbezierpath(BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(𝒑,K[i],K[i+1]))...) for i ∈ 1:N]),:stroke)

    if (points)
        sethue("black")
        setline(zoom)
        CtrlPts=[LxrPt(a[i,:],step) for i ∈ 1:size(a)[1]]
        map(p->circle(p,3*zoom,:fill), CtrlPts)
        poly(CtrlPts[:], :stroke)
    end

    finish()
    ChangeUnit(filename,"pt",unit)
    return nothing
end

function NURBSdraw(N2::NURBS2mfd;filename="NURBS2mfd.svg",up=5,down=-5,right=5,left=-5,zoom=1,mesh=(10,10),partitionnumber=(5,5),unitlength=(100,"pt"),points=true)
    step, unit=(unitlength[1]*zoom,unitlength[2])
    Drawing(step*(right-left),step*(up-down),filename)

    Luxor.origin(-step*left,step*up)
    setline(zoom)
    background("white")

    p₁,p₂=p=N2.p
    k₁,k₂=k=N2.k
    w=N2.w
    a=N2.a
    n₁,n₂=n=length.(k)-p.-1
    𝒑(u)=NURBSmapping(N2,u)

    K₁=DelDpl(k₁[1+p₁:end-p₁])
    KK=Float64[]
    for i in 1:length(K₁)-1
        for j in 0:partitionnumber[1]
            push!(KK,((partitionnumber[1]-j)*K₁[i]+j*K₁[i+1])/partitionnumber[1])
        end
    end
    K₁=DelDpl(KK)
    K₂=DelDpl(k₂[1+p₂:end-p₂])
    KK=Float64[]
    for i in 1:length(K₂)-1
        for j in 0:partitionnumber[2]
            push!(KK,((partitionnumber[2]-j)*K₂[i]+j*K₂[i+1])/partitionnumber[2])
        end
    end
    K₂=DelDpl(KK)

    K=[K₁,K₂]
    N₁,N₂=length.(K).-1
    m₁,m₂=mesh

    sethue(1,.5,.5)
    drawbezierpath(BezierPath(vcat(
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u₁->𝒑([u₁,K₂[1]]),K₁[i],K₁[i+1]))...) for i ∈ 1:N₁],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u₂->𝒑([K₁[end],u₂]),K₂[i],K₂[i+1]))...) for i ∈ 1:N₂],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u₁->𝒑([u₁,K₂[end]]),K₁[end-i+1],K₁[end-i]))...) for i ∈ 1:N₁],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u₂->𝒑([K₁[1],u₂]),K₂[end-i+1],K₂[end-i]))...) for i ∈ 1:N₂]
    )),:fill,close=true)

    sethue("red")
    for u₁ ∈ range(K₁[1],stop=K₁[end],length=m₁+1)
        drawbezierpath(BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(u₂->𝒑([u₁,u₂]),K₂[i],K₂[i+1]))...) for i ∈ 1:N₂]),:stroke)
    end
    for u₂ ∈ range(K₂[1],stop=K₂[end],length=m₂+1)
        drawbezierpath(BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(u₁->𝒑([u₁,u₂]),K₁[i],K₁[i+1]))...) for i ∈ 1:N₁]),:stroke)
    end

    if (points)
        sethue("black")
        setline(zoom)
        CtrlPts=[LxrPt(a[i,j,:],step) for i ∈ 1:size(a)[1], j ∈ 1:size(a)[2]]
        map(p->circle(p,3*zoom,:fill), CtrlPts)
        for i ∈ 1:n₁
            poly(CtrlPts[i,:], :stroke)
        end
        for j ∈ 1:n₂
            poly(CtrlPts[:,j], :stroke)
        end
    end

    finish()
    ChangeUnit(filename,"pt",unit)
    return nothing
end

function pref(N2::NURBS2mfd,p₊::Array{Int64,1})
    p,k,w,a = N2.p,N2.k,N2.w,N2.a
    p₁,p₂=p
    k₁,k₂=k
    p₁′,p₂′=p′=p+p₊
    k₊=[repeat(DelDpl(k[l]),inner=p₊[l]) for l ∈ 1:2]
    k₁′,k₂′=k′=[sort(convert(Array{Float64,1},vcat(k[l],k₊[l]))) for l ∈ 1:2]

    n₁,n₂=n=length.(k)-p.-1
    n₁′,n₂′=n′=length.(k′)-p′.-1
    A₁=hcat([BsCoef(t->Bs(i,p₁,k₁,t),p₁′,k₁′) for i in 1:n₁]...)
    A₂=hcat([BsCoef(t->Bs(i,p₂,k₂,t),p₂′,k₂′) for i in 1:n₂]...)
    w′=[sum(w[I₁,I₂]*A₁[J₁,I₁]*A₂[J₂,I₂] for I₁ in 1:n₁, I₂ in 1:n₂) for J₁ in 1:n₁′, J₂ in 1:n₂′]
    a′=[sum(a[I₁,I₂,j]*w[I₁,I₂]*A₁[J₁,I₁]*A₂[J₂,I₂]/w′[J₁,J₂] for I₁ in 1:n₁, I₂ in 1:n₂) for J₁ in 1:n₁′, J₂ in 1:n₂′, j in 1:2]
    return NURBS2mfd(p′,k′,w′,a′)
end

function href(N2::NURBS2mfd,k₊::Array{Array{Float64,1},1})
    p,k,w,a = N2.p,N2.k,N2.w,N2.a
    p₁,p₂=p
    k₁,k₂=k
    p₁′,p₂′=p′=p
    k₁′,k₂′=k′=[sort(vcat(k[l],k₊[l])) for l ∈ 1:2]

    n₁,n₂=n=length.(k)-p.-1
    n₁′,n₂′=n′=length.(k′)-p′.-1
    A₁=hcat([BsCoef(t->Bs(i,p₁,k₁,t),p₁′,k₁′) for i in 1:n₁]...)
    A₂=hcat([BsCoef(t->Bs(i,p₂,k₂,t),p₂′,k₂′) for i in 1:n₂]...)
    w′=[sum(w[I₁,I₂]*A₁[J₁,I₁]*A₂[J₂,I₂] for I₁ in 1:n₁, I₂ in 1:n₂) for J₁ in 1:n₁′, J₂ in 1:n₂′]
    a′=[sum(a[I₁,I₂,j]*w[I₁,I₂]*A₁[J₁,I₁]*A₂[J₂,I₂]/w′[J₁,J₂] for I₁ in 1:n₁, I₂ in 1:n₂) for J₁ in 1:n₁′, J₂ in 1:n₂′, j in 1:2]
    return NURBS2mfd(p′,k′,w′,a′)
end

function N(N2::NURBS2mfd,I₁,I₂,u)
    p₁,p₂=p=N2.p
    k₁,k₂=k=N2.k
    w=N2.w
    n₁,n₂=length.(k)-p.-1
    return Bs(I₁,p₁,k₁,u[1])*Bs(I₂,p₂,k₂,u[2])*w[I₁,I₂]/sum(Bs(J₁,p₁,k₁,u[1])*Bs(J₂,p₂,k₂,u[2])*w[J₁,J₂] for J₁ ∈ 1:n₁, J₂ ∈ 1:n₂)
end

function N′(N2::NURBS2mfd,I₁,I₂,i,u)
    p₁,p₂=p=N2.p
    k₁,k₂=k=N2.k
    w=N2.w
    n₁,n₂=length.(k)-p.-1
    if(i==1)
        return sum(-Bs(I₁,p₁,k₁,u[1])*Bs(I₂,p₂,k₂,u[2])*w[I₁,I₂]*Ḃs(J₁,p₁,k₁,u[1])*Bs(J₂,p₂,k₂,u[2])*w[J₁,J₂]+Ḃs(I₁,p₁,k₁,u[1])*Bs(I₂,p₂,k₂,u[2])*w[I₁,I₂]*Bs(J₁,p₁,k₁,u[1])*Bs(J₂,p₂,k₂,u[2])*w[J₁,J₂] for J₁ ∈ 1:n₁, J₂ ∈ 1:n₂)/((sum(Bs(J₁,p₁,k₁,u[1])*Bs(J₂,p₂,k₂,u[2])*w[J₁,J₂] for J₁ ∈ 1:n₁, J₂ ∈ 1:n₂))^2)
    else
        return sum(-Bs(I₁,p₁,k₁,u[1])*Bs(I₂,p₂,k₂,u[2])*w[I₁,I₂]*Bs(J₁,p₁,k₁,u[1])*Ḃs(J₂,p₂,k₂,u[2])*w[J₁,J₂]+Bs(I₁,p₁,k₁,u[1])*Ḃs(I₂,p₂,k₂,u[2])*w[I₁,I₂]*Bs(J₁,p₁,k₁,u[1])*Bs(J₂,p₂,k₂,u[2])*w[J₁,J₂] for J₁ ∈ 1:n₁, J₂ ∈ 1:n₂)/((sum(Bs(J₁,p₁,k₁,u[1])*Bs(J₂,p₂,k₂,u[2])*w[J₁,J₂] for J₁ ∈ 1:n₁, J₂ ∈ 1:n₂))^2)
    end
end

function N′(N2::NURBS2mfd,I₁,I₂,u)
    p₁,p₂=p=N2.p
    k₁,k₂=k=N2.k
    w=N2.w
    n₁,n₂=length.(k)-p.-1
    return -w[I₁,I₂]*sum(
        [Bs(I₁,p₁,k₁,u[1])*Bs(I₂,p₂,k₂,u[2])*Ḃs(J₁,p₁,k₁,u[1])*Bs(J₂,p₂,k₂,u[2])-Ḃs(I₁,p₁,k₁,u[1])*Bs(I₂,p₂,k₂,u[2])*Bs(J₁,p₁,k₁,u[1])*Bs(J₂,p₂,k₂,u[2])
         Bs(I₁,p₁,k₁,u[1])*Bs(I₂,p₂,k₂,u[2])*Bs(J₁,p₁,k₁,u[1])*Ḃs(J₂,p₂,k₂,u[2])-Bs(I₁,p₁,k₁,u[1])*Ḃs(I₂,p₂,k₂,u[2])*Bs(J₁,p₁,k₁,u[1])*Bs(J₂,p₂,k₂,u[2])]*w[J₁,J₂]
        for J₁ ∈ 1:n₁, J₂ ∈ 1:n₂)/sum(Bs(J₁,p₁,k₁,u[1])*Bs(J₂,p₂,k₂,u[2])*w[J₁,J₂] for J₁ ∈ 1:n₁, J₂ ∈ 1:n₂)^2
end

function NURBSṁapping(N2::NURBS2mfd,u)
    p₁,p₂=p=N2.p
    k₁,k₂=k=N2.k
    a=N2.a
    n₁,n₂=n=length.(k)-p.-1
    return [sum(N′(N2,I₁,I₂,j,u)*a[I₁,I₂,i] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂) for i in 1:2, j in 1:2]
end

# function NURBSjacobian(N2::NURBS2mfd,u)
#     det(NURBSṁapping(N2,u))
# end

end
