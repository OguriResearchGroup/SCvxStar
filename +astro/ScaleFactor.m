classdef ScaleFactor < handle

    properties
        
        l
        t
        mass
        v
        a
        X
        stocacc
        force
        GM

    end

    properties (Dependent)
        TU2DAYS
        DAYS2TU
    end

    methods

        function obj = ScaleFactor(l, t, mass)

            arguments
                l
                t
                mass = [];
            end

            obj.l = l;
            obj.t = t;
            obj.mass = mass;

            obj.v       = obj.l/obj.t;
            obj.a       = obj.v/obj.t;
            obj.stocacc = obj.l / obj.t^(3/2);
            obj.force   = obj.mass * obj.a;
            obj.GM      = obj.l^3 / obj.t^2;

        end

        function TU2DAYS = get.TU2DAYS(obj)
            TU2DAYS = obj.t / (24*3600);
        end

        function DAYS2TU = get.DAYS2TU(obj)
            DAYS2TU = (24*3600) / obj.t;
        end

    end
end