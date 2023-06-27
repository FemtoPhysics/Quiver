const std = @import("std");

pub fn readdlm(pth: []const u8, comptime T: type, comptime row_sz: usize, comptime col_sz: usize) anyerror![row_sz][col_sz]T {
    const file = try std.fs.cwd().openFile(pth, .{ .mode = .read_only });
    defer file.close();

    var buffer: [1024]u8 = undefined;
    const read_all_index: usize = try file.readAll(&buffer);

    var data: [row_sz][col_sz]T = undefined;
    var left: usize = 0;
    var xind: usize = 0;
    var yind: usize = 0;

    for (buffer[0..read_all_index], 0..) |char, right| {
        switch (char) {
            32, 43, 45, 46, 48...57, 69, 101 => continue,
            '\t', ',' => {
                data[0][xind] = try std.fmt.parseFloat(T, buffer[left..right]);
                left = right + 1;
                xind += 1;
            },
            '\n' => {
                data[1][yind] = try std.fmt.parseFloat(T, buffer[left..right]);
                left = right + 1;
                yind += 1;
            },
            else => unreachable,
        }
    }
    return data;
}

test "readdlm" {
    var data = try readdlm("../assets/data-2x12.tsv", f128, 2, 12);
    try std.testing.expectEqual(data[0], .{ 0.07, 0.51, 1.21, 1.40, 2.03, 2.51, 3.07, 3.43, 3.95, 4.44, 5.05, 5.38 });
    try std.testing.expectEqual(data[1], .{ 3.74, 3.96, 4.30, 4.40, 4.72, 4.96, 5.24, 5.42, 5.68, 5.92, 6.22, 6.39 });
}
