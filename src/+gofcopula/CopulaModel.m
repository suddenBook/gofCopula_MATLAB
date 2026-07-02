classdef CopulaModel
    %COPULAMODEL Immutable-style specification of a fitted or estimable copula.

    properties (SetAccess = immutable)
        Family (1,1) string
        Theta double
        DegreesOfFreedom (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 4
        EstimateTheta (1,1) logical = true
        EstimateDegreesOfFreedom (1,1) logical = true
        Dispersion (1,1) string {mustBeMember(Dispersion,["exchangeable","unstructured"])} = "exchangeable"
        Rotation (1,1) double {mustBeMember(Rotation,[0,90,180,270])} = 0
    end

    methods
        function obj = CopulaModel(family, options)
            arguments
                family {mustBeTextScalar}
                options.Theta {mustBeNumeric,mustBeReal,mustBeFinite} = 0.5
                options.DegreesOfFreedom (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBePositive} = 4
                options.EstimateTheta (1,1) logical = true
                options.EstimateDegreesOfFreedom (1,1) logical = true
                options.Dispersion {mustBeTextScalar} = "exchangeable"
                options.Rotation (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger} = 0
            end
            family = lower(string(family));
            if family == "gaussian"
                family = "normal";
            end
            mustBeMember(family, gofcopula.internal.utilities.families());
            dispersion = validatestring(options.Dispersion, ...
                {'exchangeable','unstructured'});
            mustBeMember(options.Rotation, [0,90,180,270]);
            obj.Family = family;
            obj.Theta = double(options.Theta(:).');
            obj.DegreesOfFreedom = double(options.DegreesOfFreedom);
            obj.EstimateTheta = options.EstimateTheta;
            obj.EstimateDegreesOfFreedom = options.EstimateDegreesOfFreedom;
            obj.Dispersion = string(dispersion);
            obj.Rotation = double(options.Rotation);
        end

        function obj = withFit(obj, theta, df)
            %WITHFIT Return a model carrying estimated parameters.
            arguments
                obj (1,1) gofcopula.CopulaModel
                theta {mustBeNumeric,mustBeReal,mustBeFinite}
                df (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBePositive} = obj.DegreesOfFreedom
            end
            obj = gofcopula.CopulaModel(obj.Family, Theta=theta, ...
                DegreesOfFreedom=df, EstimateTheta=false, ...
                EstimateDegreesOfFreedom=false, Dispersion=obj.Dispersion, ...
                Rotation=obj.Rotation);
        end
    end
end
