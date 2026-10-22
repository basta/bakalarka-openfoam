# backend/run_server.jl
using Pkg

# Get the directory of the run_server.jl script
script_dir = @__DIR__ # This gives the directory of run_server.jl, which is backend/

# Activate the project environment located in script_dir (i.e., backend/)
println("Activating project environment at: ", script_dir)
Pkg.activate(script_dir)

# Ensure all dependencies from Project.toml are installed and up-to-date
println("Instantiating project dependencies...")
Pkg.instantiate()

# Include the server logic
# Assuming server.jl is in backend/src/
server_file_path = joinpath(script_dir, "src", "server.jl")
println("Including server logic from: ", server_file_path)
include(server_file_path)

# Define host and port
host = "127.0.0.1"
port = 8080 # You can change this if needed

# Run the server
# The function run_makie_server should be defined in server.jl
println("Attempting to start server...")
# We call the function directly as it's now in scope after include()
# No module prefix needed unless server.jl defines one and exports the function.
# If server.jl wrapped everything in a module `MyServerModule`, you'd call `MyServerModule.run_makie_server(...)`
run_makie_server(host, port)

println("Server has been initiated. If it exited immediately, check for errors above.")
