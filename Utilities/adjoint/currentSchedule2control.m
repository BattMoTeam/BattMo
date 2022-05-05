function u = currentSchedule2control(schedule, scaling)
% Convert schedule to control vector
nc = numel(schedule.control);
[umin, umax] = deal(scaling.boxLims(:,1), scaling.boxLims(:,2));
u = cell(nc, 1);
for c = 1:nc
    ui   = vertcat(schedule.control(c).Imax);
    u{c} = (ui-umin)./(umax-umin);
end
u = vertcat(u{:});
end