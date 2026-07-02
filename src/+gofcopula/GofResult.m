classdef GofResult
    %GOFRESULT Result of one or more copula goodness-of-fit tests.

    properties (SetAccess = immutable)
        Method (1,1) string
        Copula (1,1) string
        Margins string
        MarginParameters cell
        Theta double
        DegreesOfFreedom double
        Rotation (1,1) double = 0
        Tests table
        NumericMode (1,1) string
    end

    methods
        function obj = GofResult(options)
            arguments
                options.Method {mustBeTextScalar}
                options.Copula {mustBeTextScalar}
                options.Margins string = "ranks"
                options.MarginParameters cell = {}
                options.Theta {mustBeNumeric,mustBeReal} = NaN
                options.DegreesOfFreedom {mustBeNumeric,mustBeReal} = NaN
                options.Rotation (1,1) {mustBeNumeric,mustBeReal} = 0
                options.Tests table = table(Size=[0,3], ...
                    VariableTypes=["string","double","double"], ...
                    VariableNames=["Name","PValue","Statistic"])
                options.NumericMode {mustBeTextScalar} = "corrected"
            end
            required = ["Name","PValue","Statistic"];
            if ~all(ismember(required, string(options.Tests.Properties.VariableNames)))
                error("gofcopula:GofResult:InvalidTests", ...
                    "Tests must contain Name, PValue, and Statistic variables.");
            end
            obj.Method = string(options.Method);
            obj.Copula = string(options.Copula);
            obj.Margins = options.Margins;
            obj.MarginParameters = options.MarginParameters;
            obj.Theta = double(options.Theta);
            obj.DegreesOfFreedom = double(options.DegreesOfFreedom);
            obj.Rotation = double(options.Rotation);
            obj.Tests = options.Tests;
            obj.NumericMode = string(options.NumericMode);
        end

        function disp(obj)
            for k = 1:numel(obj)
                label = obj(k).Copula;
                if obj(k).Rotation ~= 0
                    label = label + sprintf("%d", obj(k).Rotation);
                end
                fprintf('%s — %s copula (%s mode)\n', obj(k).Method, label, ...
                    obj(k).NumericMode);
                fprintf('theta: %s\n', mat2str(obj(k).Theta, 8));
                disp(obj(k).Tests);
            end
        end

        function h = plot(obj, options)
            %PLOT Plot test p-values grouped by copula family.
            arguments
                obj gofcopula.GofResult
                options.IncludeHybrid (1,1) logical = true
                options.Parent = []
            end
            if isempty(options.Parent), ax = axes(figure); else, ax = options.Parent; end
            hold(ax,"on"); labels=strings(0,1); values=[]; groups=[];
            for k=1:numel(obj)
                rows=obj(k).Tests;
                if ~options.IncludeHybrid, rows=rows(~startsWith(rows.Name,"hybrid"),:); end
                values=[values;rows.PValue]; %#ok<AGROW>
                groups=[groups;repmat(k,height(rows),1)]; %#ok<AGROW>
                labels(end+1)=obj(k).Copula; %#ok<AGROW>
            end
            h=scatter(ax,groups,values,36,"filled",MarkerFaceAlpha=0.7);
            yline(ax,0.05,"--","5% level");
            ax.XTick=1:numel(labels); ax.XTickLabel=labels; ylim(ax,[0,1]);
            ylabel(ax,"Bootstrap p-value"); xlabel(ax,"Copula"); grid(ax,"on");
            hold(ax,"off");
        end
    end
end
