import sympy as sp
from sympy.utilities.lambdify import lambdify
from math import exp, log, ceil


def input_data():
    global polynomial, start, end, x
    deg = int(input("Which degree of function?: "))
    print("Enter function")
    for i in range(deg, -1, -1):
        t = float(input(f"(x**{i})*"))
        polynomial += t * x ** i
    start = float(input("Start point: "))
    end = float(input("End point: "))
    print(f"Your func: {polynomial}, Range[{start},{end}]")


def newton_raphson(polynom, start, end):
    f = lambdify(x,polynom)
    f_dif = sp.diff(polynom,x)
    f_dif = lambdify(x, f_dif)
    counter = 0
    x_next = (start + end) / 2
    while abs(f(x_next)) > eps:
        counter += 1
        x_current = x_next
        x_next = x_current - f(x_current)/f_dif(x_current)
    return [x_next, counter]


def bisection_method(polynom, start, end):
    f = lambdify(x, polynom)
    max_iterations = log(eps/(end-start), exp(1))
    max_iterations /= -log(2, exp(1))
    max_iterations = ceil(max_iterations)
    x_l = start
    x_r = end
    counter = 0
    x_c = (x_l + x_r) / 2
    while abs(x_r - x_l) > eps:
        counter += 1
        x_c = (x_l + x_r) / 2
        if(f(x_c) * f(x_r)) < 0:
            x_l = x_c
        elif(f(x_c) * f(x_r)) > 0:
            x_r = x_c
        if counter > max_iterations:
            print("cannot resolve")
            return None
    return [x_c, counter]


def secant_method(polynom, start, end):
    f = lambdify(x, polynom)
    x_current = start
    x_prev = end
    counter = 0
    while abs(x_current - x_prev) > eps:
        counter += 1
        x_next = (x_prev*f(x_current) - x_current*f(x_prev)) / (f(x_current) - f(x_prev))
        x_prev = x_current
        x_current = x_next
    return [x_current, counter]


def main():
    global polynomial, start, end, x
    x_l = start
    print("""Choose preferred method:
    1 - Bisection method
    2 - Newton-Raphson method
    3 - Secant method""")
    choice = int(input())
    method = None
    if choice == 1:
        method = bisection_method
    elif choice == 2:
        method = newton_raphson
    elif choice == 3:
        method = secant_method
    else:
        print("Error")
    func = sp.lambdify(x, polynomial)
    p_dif = sp.diff(polynomial, x)
    f_dif = sp.lambdify(x, p_dif)
    solution = []
    while x_l < end:
        if func(x_l)*func(x_l + 0.1) < 0:
            temp = method(polynomial, x_l, (x_l + 0.1))
            if temp is not None:
                solution.append(temp)
        elif f_dif(x_l)*f_dif(x_l+0.1) < 0:
            temp = method(p_dif, x_l, (x_l + 0.1))
            if temp is not None:
                if abs(func(temp[0])) < eps:
                    solution.append(temp)
        x_l += 0.1
    for i in solution:
        if i is not None:
            print("%.5f" % i[0], end=", ")
            print(f"found in {i[1]} attempts")


x = sp.symbols('x')
eps = 10**-10
polynomial = 0
start = 0
end = 0
input_data()
main()
