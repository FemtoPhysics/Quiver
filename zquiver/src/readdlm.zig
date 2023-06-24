const std = @import("std");
const fs = std.fs;
const print = std.debug.print;

pub fn main() anyerror!void {
    const file = try fs.cwd().openFile(
        "../dev/linear-regression.dat",
        .{ .mode = .read_only }, // fs.File.OpenFlags{ .mode = .read_only }
    );
    defer file.close();

    var buffer: [1024]u8 = undefined;
    const read_all_index: usize = try file.readAll(&buffer);

    var data: [2][128]f64 = undefined;
    var left: usize = 0;
    var xind: usize = 0;
    var yind: usize = 0;

    for (buffer[0..read_all_index], 0..) |char, right| {
        switch (char) {
            32, 43, 45, 46, 48...57, 69, 101 => continue,
            '\t', ',' => {
                data[0][xind] = try std.fmt.parseFloat(f64, buffer[left..right]);
                left = right + 1;
                xind += 1;
            },
            '\n' => {
                data[1][yind] = try std.fmt.parseFloat(f64, buffer[left..right]);
                left = right + 1;
                yind += 1;
            },
            else => unreachable,
        }
    }

    print("x = {any}\n", .{data[0][0..xind]});
    print("y = {any}\n", .{data[1][0..yind]});
}
