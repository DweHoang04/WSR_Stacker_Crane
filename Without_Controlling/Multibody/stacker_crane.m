function [xdd,ydd,th1dd,th2dd,th3dd,th4dd,w] = stacker_crane(Parameters, Forces, y, th1, th2, th3, th4, yd, th1d, th2d, th3d, th4d)
    % Parameters
    L = Parameters(1);
    M = Parameters(2);
    m1 = Parameters(3);
    m2 = Parameters(4);
    m3 = Parameters(5);
    g = Parameters(6);
    k = Parameters(7);

    % Define a small space step for boundary conditions
    epsilon = 1e-6;
    
    % Input forces
    F1   = Forces(1,1);
    F2   = Forces(2,1);
    
    % Trigonometries
    s1 = sin(th1);   c1 = cos(th1);
    s2 = sin(th2);   c2 = cos(th2);
    s3 = sin(th3);   c3 = cos(th3);
    s4 = sin(th4);   c4 = cos(th4);
    
    s12 = sin(th1 - th2); c12 = cos(th1 - th2);
    s13 = sin(th1 - th3); c13 = cos(th1 - th3);
    s14 = sin(th1 - th4); c14 = cos(th1 - th4);
    s23 = sin(th2 - th3); c23 = cos(th2 - th3);
    s24 = sin(th2 - th4); c24 = cos(th2 - th4);
    s34 = sin(th3 - th4); c34 = cos(th3 - th4);
    
    % ----- STACKER CRANE DYNAMICS -----
    if y < 0
        % If y is below the mast, calculate as if it's at the very bottom (first pendulum)
        y_calc = epsilon;

        % ----- First pendulum dynamics -----
        Phi = [
            (F1 + (7*m1/32*L + m3/4*L + m2*y_calc)*s1*th1d^2 ...
                 + (5*m1/32 + m3/4)*L*s2*th2d^2);
            (F2 + m2*y_calc*th1d^2 + m2*g*c1);
           (-(5*m1/128 + m3/16)*L^2*s12*th2d^2 ...
            -(3*m1/128 + m3/16)*L^2*s13*th3d^2 ...
            -(m1/128 + m3/16)*L^2*s14*th4d^2 ...
            -2*m2*y_calc*yd*th1d - k*th1 + m2*g*s1*y_calc ...
            + (7*m1/32 + m3/4)*g*L*s1);
            ((5*m1/128 + m3/16)*L^2*s12*th1d^2 ...
            -(3*m1/128 + m3/16)*L^2*s23*th3d^2 ...
            -(m1/128 + m3/16)*L^2*s24*th4d^2 ...
            - k*th2 + (5*m1/32 + m3/4)*g*L*s2);
            ((3*m1/128 + m3/16)*L^2*(s13*th1d^2 + s23*th2d^2) ...
            -(m1/128 + m3/16)*L^2*s34*th4d^2 ...
            - k*th3 + (3*m1/32 + m3/4)*g*L*s3);
            ((m1/128 + m3/16)*L^2*(s14*th1d^2 + s24*th2d^2 + s34*th3d^2) ...
            - k*th4 + (m1/32 + m3/4)*g*L*s4)
        ];
        Mmat = [
            (M + m1 + m2 + m3), m2*s1, (7*m1/32*L + m3/4*L + m2*y_calc), (5*m1/32 + m3/4)*L*c2, (3*m1/32 + m3/4)*L*c3, (m1/32 + m3/4)*L*c4;
            m2*s1,        m2,          0,                           0,                           0,                           0;
            (m2*y_calc + 7*m1/32*L + m3/4*L)*c1, 0, (m2*y_calc^2 + 13*m1/256*L^2 + m3/16*L^2), (5*m1/128 + m3/16)*L^2*c12, (3*m1/128 + m3/16)*L^2*c13, (m1/128 + m3/16)*L^2*c14;
            (5*m1/32 + m3/4)*L*c2, 0, (5*m1/128 + m3/16)*L^2*c12, (9*m1/256 + m3/16)*L^2, (3*m1/128 + m3/16)*L^2*c23, (m1/128 + m3/16)*L^2*c24;
            (3*m1/32 + m3/4)*L*c3, 0, (3*m1/128 + m3/16)*L^2*c13, (3*m1/128 + m3/16)*L^2*c23, (5*m1/256 + m3/16)*L^2, (m1/128 + m3/16)*L^2*c34;
            (m1/32 + m3/4)*L*c4, 0, (m1/128 + m3/16)*L^2*c14, (m1/128 + m3/16)*L^2*c24, (m1/128 + m3/16)*L^2*c34, (m1/256 + m3/16)*L^2
        ];

        % Calculating outputs
        qdd = Mmat \ Phi;

        % Tip displacement
        w = y_calc*s1;

    elseif y > L
        
        % If y is above the mast, calculate as if it's at the very top (fourth pendulum)
        y_calc = L - epsilon;

        % ----- Fourth pendulum dynamics -----
        k1 = M + m1 + m2 + m3;
        k2 = (7*m1/32 + (m2 + m3)/4)*L;
        k3 = (5*m1/32 + (m2 + m3)/4)*L;
        k4 = (3*m1/32 + (m2 + m3)/4)*L;
        k5 = m1*L/32 - 3*m2*L/4 + m3*L/4 + m2*y_calc;
        k6 = (5*m1/128 + (m2 + m3)/16)*L^2;
        k7 = (3*m1/128 + (m2 + m3)/16)*L^2;
        k8 = 13*m1*L^2/128 - 3*m2*L^2/16 + m3*L^2/16 + m2*y_calc*L/4;
        k9 = m1*L^2/128 - 3*m2*L^2/16 + m3*L^2/16 + m2*y_calc*L/4;
        k10 = m1*L^2/128 - (3*m2 + m3)*L^2/16 + m2*y_calc*L/4;
        Phi = [
            (F1 + k2*s1*th1d^2 + k3*s2*th2d^2 + k4*s3*th3d^2 + k5*s4*th4d^2 - 2*m2*c4*yd*th4d);
            (F2 + m2*L/4*(c14*th1d^2 + c24*th2d^2 + c34*th3d^2) - (3*m2*L/4 - m2*y_calc)*th4d^2 + m2*g*c4);
            (-k6*s12*th2d^2 - k7*s13*th3d^2 - k8*s14*th4d^2 - m2*L/2*c14*yd*th4d - k*th1 + (7*m1/32 + (m2 + m3)/4)*g*L*s1);
            ((5*m1/128 + (m2 + m3)/16)*s12*th1d^2 - k7*s23*th3d^2 - m2*L*c24*yd*th4d + k*th1 - k9*s24*th4d^2 - k*th2 + k3*g*s2);
            (k7*(s13*th1d^2 + s23*th2d^2) - k10*s34*th4d^2 - m2*L/2*c34*yd*th4d - k*th3 + k4*g*s3);
            (-(2*m2*y_calc - 3*m2*L/2)*yd*th4d + k9*(s14*th1d^2 + s24*th2d^2 + s34*th3d^2) - k*th4 + (m2*y_calc + m1*L/32 - 3*m2*L/4 + m3*L/4)*g*s4)
        ];
        Mmat = [
            k1, m2*s4, k2*c1, k3*c2, k4*c3, k5*c4;
            m2*s4, m2, -m2*L*s14/4, -m2*L*s24/4, -m2*L*s34/4, 0;
            k2*c1, -m2*L*s14/4, (13*m1/256 + m2/16 + m3/16)*L^2, k6*c12, k7*c13, k8*c14;
            k3*c2, -m2*L*s24/4, k6*c12, (9*m1/256 + (m2 + m3)/16)*L^2, k7*c23, k9*c24;
            k4*c3, -m2*L*s34/4, k7*c13, k7*c23, (5*m1/256 + (m2 + m3)/16)*L^2, k10*c34;
            (m2*y_calc + m1*L/32 - 3*m2*L/4 + m3*L/4)*c4, 0, k9*c14, k9*c24, k9*c34, (m2*y_calc^2 + m1*L^2/256 + (9*m2 + m3)*L^2/16 - 3*m2*y_calc*L/2)
        ];
        
        % Calculating outputs
        qdd = Mmat \ Phi;
        
        % Tip displacement
        w = (L/4)*(s1 + s2 + s3) + (y_calc - 3*L/4)*s4;

    else
        % --- ORIGINAL LOGIC FOR 0 <= y <= L ---
        if y >= 0 && y < (L/4)
            % ----- First pendulum -----
            Phi = [
                (F1 + (7*m1/32*L + m3/4*L + m2*y)*s1*th1d^2 ...
                     + (5*m1/32 + m3/4)*L*s2*th2d^2);
                (F2 + m2*y*th1d^2 + m2*g*c1);
               (-(5*m1/128 + m3/16)*L^2*s12*th2d^2 ...
                -(3*m1/128 + m3/16)*L^2*s13*th3d^2 ...
                -(m1/128 + m3/16)*L^2*s14*th4d^2 ...
                -2*m2*y*yd*th1d - k*th1 + m2*g*s1*y ...
                + (7*m1/32 + m3/4)*g*L*s1);
                ((5*m1/128 + m3/16)*L^2*s12*th1d^2 ...
                -(3*m1/128 + m3/16)*L^2*s23*th3d^2 ...
                -(m1/128 + m3/16)*L^2*s24*th4d^2 ...
                - k*th2 + (5*m1/32 + m3/4)*g*L*s2);
                ((3*m1/128 + m3/16)*L^2*(s13*th1d^2 + s23*th2d^2) ...
                -(m1/128 + m3/16)*L^2*s34*th4d^2 ...
                - k*th3 + (3*m1/32 + m3/4)*g*L*s3);
                ((m1/128 + m3/16)*L^2*(s14*th1d^2 + s24*th2d^2 + s34*th3d^2) ...
                - k*th4 + (m1/32 + m3/4)*g*L*s4)
            ];
            Mmat = [
                (M + m1 + m2 + m3), m2*s1, (7*m1/32*L + m3/4*L + m2*y), (5*m1/32 + m3/4)*L*c2, (3*m1/32 + m3/4)*L*c3, (m1/32 + m3/4)*L*c4;
                m2*s1,        m2,          0,                           0,                           0,                           0;
                (m2*y + 7*m1/32*L + m3/4*L)*c1, 0, (m2*y^2 + 13*m1/256*L^2 + m3/16*L^2), (5*m1/128 + m3/16)*L^2*c12, (3*m1/128 + m3/16)*L^2*c13, (m1/128 + m3/16)*L^2*c14;
                (5*m1/32 + m3/4)*L*c2, 0, (5*m1/128 + m3/16)*L^2*c12, (9*m1/256 + m3/16)*L^2, (3*m1/128 + m3/16)*L^2*c23, (m1/128 + m3/16)*L^2*c24;
                (3*m1/32 + m3/4)*L*c3, 0, (3*m1/128 + m3/16)*L^2*c13, (3*m1/128 + m3/16)*L^2*c23, (5*m1/256 + m3/16)*L^2, (m1/128 + m3/16)*L^2*c34;
                (m1/32 + m3/4)*L*c4, 0, (m1/128 + m3/16)*L^2*c14, (m1/128 + m3/16)*L^2*c24, (m1/128 + m3/16)*L^2*c34, (m1/256 + m3/16)*L^2
            ];
            
            % Calculating outputs
            qdd = Mmat \ Phi;

            % Tip displacement
            w = y*s1;
    
        elseif y >= (L/4) && y < (L/2)
            % ----- Second pendulum -----
            k1 = M + m1 + m2 + m3;
            k2 = (7*m1/32 + (m2 + m3)/4)*L;
            k3 = 5*L*m1/32 - (m2 - m3)*L/4 + m2*y;
            k4 = (3*m1/32 + m3/4)*L;
            k5 = (m1/32 + m3/4)*L;
            k6 = (5*m1*L^2/128 + (m2*L^2 + m3*L^2)/16 + m2*L*y/4);
            k7 = (13*m1/128 + m3/16)*L^2;
            k8 = 5*m1*L^2/128 - (m2 - m3)*L^2/16 + m2*L*y/4;
            k9 = (3*m1/128 + m3/16)*L^2;
            k10 = (m1/128 + m3/16)*L^2;
            Mmat = [
                k1,             m2*s2,                        k2*c1,                      k3*c2,                 k4*c3,                     k5*c4;
                m2*s2,         m2,                           -m2*L/4*s12,                 0,                     0,                         0;
                k2*c1,        -m2*L/4*s12,   (13*m1/256+(m2+m3)/16)*L^2,   k6*c12,     k7*c13,                 k7*c14;
                k3*c2,         k8*c12,         0,           (9*m1*L^2/256 + m2*L^2/16 + m3*L^2/16 + m2*y^2 - m2*L*y/2), k9*c23, k10*c24;
                k4*c3,         0,             k9*c13,        k9*c23,       (5*m1/256 + m3/16)*L^2,              k10*c34;
                k5*c4,         0,             k10*c14,       k10*c24,      k10*c34,            (m1/256 + m3/16)*L^2
            ];
            Phi = [
                (F1 + k2*s1*th1d^2 + k3*s2*th2d - 2*m2*c2*yd*th2d + k4*s3*th3d^2 + k5*s4*th4d^2);
                (F2 + m2*L/4*c12*th1d^2 + m2*g*c2 + (-m2*L/4 + m2*y)*th2d^2);
                (-k6*s12*th2d^2 - m2*L/2*c12*yd*th2d - k7*(s13*th3d^2 + s14*th4d^2) - k*th1 + (7*m1/32+(m2+m3)/4)*g*L*s1);
                (k8*s12*th1d^2 - (2*m2*y - m2*L/2)*yd*th2d - k9*s23*th3d^2 - k10*s24*th4d^2 - k*th2 + (m2*y+5*m1*L/32-(m2-m3)*L/4)*g*s2);
                ((3*m1/128+m3/16)*L^2*(s13*th1d^2+s23*th2d^2) - (m1/128+m3/16)*L^2*s34*th4d^2 - k*th3 + (3*m1/32+m3/4)*g*L*s3);
                ((m1/128+m3/16)*L^2*(s14*th1d^2+s24*th2d^2+s34*th3d^2) - k*th4 + (m1/32+m3/4)*g*L*s4)
            ];
            
            % Calculating outputs
            qdd = Mmat \ Phi;

            % Tip displacement
            w = (L/4)*s1 + (y - L/4)*s2;
    
        elseif y >= (L/2) && y < (3*L/4)

            % ----- Third pendulum -----
            k1  = M + m1 + m2 + m3;
            k2  = (7*m1/32 + (m2 + m3)/4)*L;
            k3  = (5*m1/32 + (m2 + m3)/4)*L;
            k4  = 3*m1*L/32 - m2*L/2 + m3*L/4 + m2*y;
            k5  = (3*m1/32 + m3/4)*L;
            k6  = m2*L/4;
            k7  = (5*m1/128 + (m2 + m3)/16)*L^2;
            k8  = 3*m1*L^2/128 - m2*L^2/8 + m3*L^2/16 + m2*L*y/4;
            k9  = (m1/128 + m3/16)*L^2;
            k10 = 5*m1/128 + m2/16 + m3/16;
            k11 = 3*m1*L^2/128 - m2*L^2/8 + m3*L^2/16 + m2*L*y/4;
            k12 = 3*m1*L^2/128 - m2*L^2/8 + m3*L^2/16 + m3*L*y/4;
            Phi = [
                (F1 + k2*s1*th1d^2 + k3*s2*th2d^2 + k4*s3*th3d^2 - 2*m2*c3*yd*th3d + k5*s4*th4d^2);
                (F2 + k6*(c13*th1d^2 + c23*th2d^2) - (m2*L/2 - m2*y)*th3d^2 + m2*g*c3);
                (-k7*s12*th2d^2 - k12*s13*th3d^2 - k9*s14*th4d^2 - k*th1 - m2*L/2*c13*yd*th3d + (7*m1/32 + (m2 + m3)/4)*g*L*s1);
                (k10*s12*th1d^2 - k*th2 - k11*s23*th3d^2 - (m1/128 + m3/16)*L^2*s24*th4d^2 - m2*L/2*c23*yd*th3d + k3*g*s2);
                (k8*(s13*th1d^2 + s23*th2d^2) - k9*s34*th4d^2 - k*th3 - (2*m2*y - m2*L)*yd*th3d + (m2*y + 3*m1*L/32 - m2*L/2 + m3*L/4)*g*s3);
                ((m1/128 + m3/16)*L^2*(s14*th1d^2 + s24*th2d^2 + s34*th3d^2) - k*th4 + (m1/32 + m3/4)*g*L*s4)
            ];
            Mmat = [
                k1, m2*s3, k2*c1, k3*c2, k4*c3, k5*c4;
                m2*s3, m2, -m2*L*s13/4, -m2*L*s23/4, 0, 0;
                k2*c1, -m2*L*s13/4, (13*m1/256 + m2/16 + m3/16)*L^2, k7*c12, k12*c13, k9*c14;
                k3*c2, -m2*L*s23/4, k10*c12, (9*m1/256 + (m2 + m3)/16)*L^2, k11*c23, k9*c24;
                (m2*y + 3*m1*L/32 - m2*L/2 + m3*L/4)*c3, 0, k8*c13, k8*c23, (m2*y^2 + 5*m1*L^2/256 + m2*L^2/4 + m3*L^2/16 - m2*L*y), k9*c34;
                (m1/32 + m3/4)*L*c4, 0, (m1/128 + m3/16)*L^2*c14, (m1/128 + m3/16)*L^2*c24, (m1/128 + m3/16)*L^2*c34, (m1/256 + m3/16)*L^2
            ];

            % Calculating outputs
            qdd = Mmat \ Phi;

            % Tip displacement
            w = (L/4)*(s1 + s2) + (y - L/2)*s3;
    
        else

            % ----- Fourth pendulum -----
            k1 = M + m1 + m2 + m3;
            k2 = (7*m1/32 + (m2 + m3)/4)*L;
            k3 = (5*m1/32 + (m2 + m3)/4)*L;
            k4 = (3*m1/32 + (m2 + m3)/4)*L;
            k5 = m1*L/32 - 3*m2*L/4 + m3*L/4 + m2*y;
            k6 = (5*m1/128 + (m2 + m3)/16)*L^2;
            k7 = (3*m1/128 + (m2 + m3)/16)*L^2;
            k8 = 13*m1*L^2/128 - 3*m2*L^2/16 + m3*L^2/16 + m2*y*L/4;
            k9 = m1*L^2/128 - 3*m2*L^2/16 + m3*L^2/16 + m2*y*L/4;
            k10 = m1*L^2/128 - (3*m2 + m3)*L^2/16 + m2*y*L/4;
            Phi = [
                (F1 + k2*s1*th1d^2 + k3*s2*th2d^2 + k4*s3*th3d^2 + k5*s4*th4d^2 - 2*m2*c4*yd*th4d);
                (F2 + m2*L/4*(c14*th1d^2 + c24*th2d^2 + c34*th3d^2) - (3*m2*L/4 - m2*y)*th4d^2 + m2*g*c4);
                (-k6*s12*th2d^2 - k7*s13*th3d^2 - k8*s14*th4d^2 - m2*L/2*c14*yd*th4d - k*th1 + (7*m1/32 + (m2 + m3)/4)*g*L*s1);
                ((5*m1/128 + (m2 + m3)/16)*s12*th1d^2 - k7*s23*th3d^2 - m2*L*c24*yd*th4d + k*th1 - k9*s24*th4d^2 - k*th2 + k3*g*s2);
                (k7*(s13*th1d^2 + s23*th2d^2) - k10*s34*th4d^2 - m2*L/2*c34*yd*th4d - k*th3 + k4*g*s3);
                (-(2*m2*y - 3*m2*L/2)*yd*th4d + k9*(s14*th1d^2 + s24*th2d^2 + s34*th3d^2) - k*th4 + (m2*y + m1*L/32 - 3*m2*L/4 + m3*L/4)*g*s4)
            ];
            Mmat = [
                k1, m2*s4, k2*c1, k3*c2, k4*c3, k5*c4;
                m2*s4, m2, -m2*L*s14/4, -m2*L*s24/4, -m2*L*s34/4, 0;
                k2*c1, -m2*L*s14/4, (13*m1/256 + m2/16 + m3/16)*L^2, k6*c12, k7*c13, k8*c14;
                k3*c2, -m2*L*s24/4, k6*c12, (9*m1/256 + (m2 + m3)/16)*L^2, k7*c23, k9*c24;
                k4*c3, -m2*L*s34/4, k7*c13, k7*c23, (5*m1/256 + (m2 + m3)/16)*L^2, k10*c34;
                (m2*y + m1*L/32 - 3*m2*L/4 + m3*L/4)*c4, 0, k9*c14, k9*c24, k9*c34, (m2*y^2 + m1*L^2/256 + (9*m2 + m3)*L^2/16 - 3*m2*L*y/2)
            ];

            % Calculating outputs
            qdd = Mmat \ Phi;

            % Tip displacement
            w = (L/4)*(s1 + s2 + s3) + (y - 3*L/4)*s4;
        end
    end
    
    % Assign final outputs
    xdd   = qdd(1);
    ydd   = qdd(2);
    th1dd = qdd(3);
    th2dd = qdd(4);
    th3dd = qdd(5);
    th4dd = qdd(6);

end