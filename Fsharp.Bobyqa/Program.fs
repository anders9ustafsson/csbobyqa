open System
open Cureos.Numerics

let hs05 (n : int) (x : float[]) : float =
    Math.Sin(x.[0] + x.[1]) + Math.Pow(x.[0] - x.[1], 2.0) - 1.5 * x.[0] + 2.5 * x.[1] + 1.0

let x = [| 0.0; 0.0 |]
let status = Bobyqa.FindMinimum(new Func<int, float[], float>(hs05), 2, x, [| -1.5; -3.0 |], [| 4.0; 3.0 |])

printfn "Optimization status: %O" status
printfn "x* = %A" x

printfn "\nPress <ENTER> to exit..."
Console.ReadLine() |> ignore
