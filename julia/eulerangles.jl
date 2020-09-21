
function euler(ϕ, θ, ψ)
    Rϕ = k_rotation(ϕ)
    Rθ = i_rotation(θ)
    Rψ = k_rotation(ψ)
    
    R = Rψ * Rθ * Rϕ
end

function i_rotation(α)
    C = cos(α)
    S = sin(α)
    R = [1 0 0; 0 C S; 0 -S C]
end

function j_rotation(α)
    C = cos(α)
    S = sin(α)
    R = [C 0 S; 0 1 0; -S 0 C]
end

function k_rotation(α)
    C = cos(α)
    S = sin(α)
    R = [C S 0; -S C 0; 0 0 1]
end

##   T E S T     F U N C T I O N S

function det2d(a)
    a[1,1] * a[2,2] - a[1,2] * a[2,1]
end

function eulertest(R, failed, passed, ϵ)
    # Each row must have unit norm
    for r in 1:3
        row = R[r, :]
        row_norm2  = sum(row .* row)
        if abs(row_norm2 -1) > ϵ
            failed[1] += 1
        else
            passed[1] += 1
        end
    end
    # Each column must have unit norm
    for c in 1:3
        col = R[:, c]
        col_norm2 = sum(col .* col)
        if abs(col_norm2 -1) > ϵ
            failed[2] += 1
        else
            passed[2] += 1
        end
    end
    # R determinant must be ± 1
    #   expanded along the 1st row
    sub1 = [R[2,2] R[2,3]; R[3,2] R[3,3]]   
    sub2 = [R[2,1] R[2,3]; R[3,1] R[3,3]]
    sub3 = [R[2,1] R[2,2]; R[3,1] R[3,2]]
    det  = R[1,1] * det2d(sub1) - R[1,2] * det2d(sub2) + R[1,3] * det2d(sub3)

    det_norm2 = det*det
    if abs(det_norm2 -1) > ϵ
        failed[3] += 1
    else
        passed[3] += 1
    end

    return det
end


##   TEST DATA 
n_ϕ = 11; d_ϕ = 2π / n_ϕ 
n_θ = 11; d_θ = 2π / n_θ 
n_ψ = 11; d_ψ = 2π / n_ψ 

ϕ_data = [k * d_ϕ for k in 0:n_ϕ]
θ_data = [k * d_θ for k in 0:n_θ]
ψ_data = [k * d_ψ for k in 0:n_ψ]

##   T E S T     E X E C U T I O N 

failed = zeros(Int, 3)           # row, column and determinant failed, in that order
passed = zeros(Int, 3)           # row, column and determinant passed, in that order
ϵ = eps(Float32)                 # Assuming Float64 being the machine default  

right = zeros(Int,1)             # right-handed among the rotated systems   (***)
left  = zeros(Int,1)             #  left-handed among the rotated systems   (***)
for ψ in ψ_data
    for θ in θ_data
        for ϕ in ϕ_data
            R = euler(ϕ, θ, ψ)
            det = eulertest(R, failed, passed, ϵ)
            if det > 0
                right[1] += 1
            else
                left[1] += 1
            end        
        end
    end
end

println("rows        test $(failed[1]) failed, ($(passed[1]) passed)")
println("columns     test $(failed[2]) failed, ($(passed[2]) passed)")
println("determinant test $(failed[3]) failed, ($(passed[3]) passed)")
println()
println("Found $(right[1]) right-handed among the rotated systems")
println("Found $(left[1])  left-handed among the rotated systems")

# (***)
# Note. Outside of a function, scalars are not visible within a for loop,
# so by defining the variables 'right' and 'left' as scalars at lines 88-89
# we would get: 
# error "UndefVarError: right not defined" at line 96  
# error "UndefVarError: left not defined" at line 98 
