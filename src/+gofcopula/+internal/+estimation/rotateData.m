function u = rotateData(u, rotation)
%ROTATEDATA Rotate bivariate pseudo-observations as in gofCopula 0.4-3.
%   Matches R's .rotateCopula exactly (original R internal_rotateCopula.R):
%     90  -> (1-u2, u1)
%     180 -> (1-u1, 1-u2)
%     270 -> (u2, 1-u1)
arguments
    u {mustBeFloat,mustBeReal,mustBeFinite,mustBeMatrix}
    rotation (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger}
end
mustBeMember(rotation,[0,90,180,270]);
if rotation ~= 0 && size(u,2) ~= 2
    error("gofcopula:Rotation:Dimension", "Rotations require bivariate data.");
end
switch rotation
    case 90,  u = [1-u(:,2), u(:,1)];
    case 180, u = 1-u;
    case 270, u = [u(:,2), 1-u(:,1)];
end
end
