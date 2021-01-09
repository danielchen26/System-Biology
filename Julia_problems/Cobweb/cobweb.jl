using Plots;gr()
r = 0.7
f(x) = r*x.*(3-x)
f(x) = x*(1 + 2.5(1-x))
function check_intercept()
    plot(0:0.1:2, f.(0:0.1:2))
    plot!(x->x,0:2)
end 
check_intercept()


function cobweb(f, a, b, N, x0)
    rr = range(a, b,length = 100) # or LinRange(a,b,N)
    plot(rr,f.(rr))
    # Plot orbit starting at x0
    plt1 = scatter!([x0],[x0],c=:red, markerstrokecolor=:red)
    x = zeros(N);x[1] = x0; [x[i+1]=f(x[i]) for i in 1:N-1];
    @show x
    for i=1:N-1
        plot!([x[i],x[i]],[x[i],x[i+1]],c=:red,legend = false)
        plot!([x[i],x[i+1]],[x[i+1],x[i+1]],c=:green,legend = false)
    end 
    # xlims!(a,b)
    # ylims!(a,b)
    plot!(x -> x, 0, 2)
    title!("Cobweb")
    return x, display(plt1)
end

##
a = 0
b=2
N =30
x0=0.05
x, _ = cobweb(f, a, b, N, x0)
plot(x)
##

