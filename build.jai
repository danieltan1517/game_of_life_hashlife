#import "Basic";
#import "Compiler";

build :: () {
    w := compiler_create_workspace();
    defer set_build_options_dc(.{do_output=false});
    if !w {
        print("Workspace creation failed.\n");
        return;
    }

    target_options := get_build_options(w);
    args := target_options.compile_time_command_line;
    gui := false;
    cmd := false;
    test := false;
    release := true;

    for arg: args {
        if arg == {

        case "gui";
            gui = true;

        case "cmd";
            cmd = true;

        case "test";
            test = true;

        case "release";
            release = true;

        case "r";
            release = true;

        case "d";
            release = false;

        case "debug";
            release = false;

        case "help";
            print("Command-line options:\n");
            print("gui         Builds the GUI Interface\n");
            print("cmd         Builds the Command-Line Interface\n");
            print("test        Builds the robustness/soak test harness (gol-test)\n");
            print("Run 'jai build.jai - gui cmd' to build both the gui and commandline application\n");
            return;
        case;

        }
    }

    if cmd == true {
        if release {
            build_release("gol-console", "hashlife_console.jai");
        } else {
            build_debug("gol-console", "hashlife_console.jai");
        }
    }
    if gui == true {
        if release {
            build_release("gol", "hashlife_ui.jai");
        } else {
            build_debug("gol", "hashlife_ui.jai");
        }
    }

    if test == true {
        // The harness is sensitive to optimization-related miscompiles, so a
        // debug build is the most useful default for catching regressions.
        // Use 'r'/'release' alongside 'test' to soak-test the optimized build.
        if release {
            build_release("gol-test", "hashlife_test.jai");
        } else {
            build_debug("gol-test", "hashlife_test.jai");
        }
    }

    if !cmd && !gui && !test {
        print("Command-line options:\n");
        print("gui         Builds the GUI Interface\n");
        print("cmd         Builds the Command-Line Interface\n");
        print("test        Builds the robustness/soak test harness (gol-test)\n");
        print("Run 'jai build.jai - gui cmd' to build both the gui and commandline application\n");
    }
}

build_debug :: (executable_name: string, build_files: ..string) -> Workspace {
    w := compiler_create_workspace(executable_name);
    options := get_build_options(w);
    options.output_executable_name = executable_name;
    set_optimization(*options, .DEBUG);
    options.backend = .X64;
    options.stack_trace = true;
    set_build_options(options, w);
    print("Creating Debug Build for: %\n", executable_name);
    for it : build_files {
        add_build_file(it, w);
    }
    return w;
}

build_release :: (executable_name: string, build_files: ..string) -> Workspace {
    w := compiler_create_workspace(executable_name);
    options := get_build_options(w);
    options.output_executable_name = executable_name;
    set_optimization(*options, .OPTIMIZED);
    options.llvm_options.enable_split_modules = false;
    options.stack_trace = false;
    set_build_options(options, w);
    print("Creating Release Build for: %\n", executable_name);
    for it : build_files {
        add_build_file(it, w);
    }
  return w;
}

#run build();
