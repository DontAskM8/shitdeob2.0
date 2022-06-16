const fs = require("fs")
const beautify = require("js-beautify")
const parser = require("@babel/parser")
const traverse = require("@babel/traverse").default
const generator = require("@babel/generator").default
const t = require("@babel/types")

var original = fs.readFileSync("./raw.js", "utf8")

var globalVars = () => Object.keys(global).filter(x => !['global', 'clearInterval', 'clearTimeout', 'setInterval', 'setTimeout', 'queueMicrotask', 'performance', 'clearImmediate', 'setImmediate'].includes(x))
function isGlobalVars(va){
	return globalVars().includes(va) && typeof global[va] != "function"
}

//Main function that triggers the entire universe
var executingFunction = /(?<=(return ))([A-z]|[0-9])+\.call\(this,([A-z]|[0-9])+\)/g
var [callingFunc, callingVal] = original.match(executingFunction)[0].replace(/call\(this,|\)/g, "").split(".")
let switchSub = []
let objName = ""
let windowName = ""

console.log(callingFunc, callingVal)

var _testCode = parser.parse(original)
traverse(_testCode, {
	//Eval all functions make it global
	"FunctionDeclaration": function(path){
		try {
			if(t.isFunctionDeclaration(path)) {
				eval(`global["${path.node.id.name}"] = ${original.slice(path.node.start, path.node.end)}`)
				if(path.node.body.body.length == 1){
					if(path.node.body.body[0].expression.type == "AssignmentExpression"){
						// console.log(path.node.body.body[0].expression)
						eval("global." + original.slice(path.node.body.body[0].expression.start, path.node.body.body[0].expression.end))
					}
				}
			}
			else{
			}
		}catch(e){
		}
	}
})
traverse(_testCode, {
	"VariableDeclaration": function(path){
		for(let declaration of path.node.declarations){
			if(objName == "" && t.isObjectExpression(declaration.init)){
				objName = declaration.id.name
				eval(`${objName} = {}`)
			}
			if(windowName == "" && t.isIdentifier(declaration.init) && declaration.init.name == "window"){
				windowName = declaration.id.name
				eval(`${windowName} = global`)
			}
			if(!t.isFunctionExpression(declaration.init)) continue;
			try {
				eval(`global.${declaration.id.name} = ${original.slice(declaration.init.start, declaration.init.end)}`)
				if(declaration.init.id){
					// var xx = function yy(){} -> yy = xx
					global[declaration.init.id.name] = eval(declaration.id.name)
				}
				
				if(declaration.init.body.body.length == 1){
					if(declaration.init.body.body[0].expression.type == "AssignmentExpression"){
						eval("global." + original.slice(declaration.init.body.body[0].expression.start, declaration.init.body.body[0].expression.end))
					}
				}
			}
			catch(e){
				// console.log(e)
			}
		}
	},
	// Get all the switch subjects so we dont accidentally remove the assignment
	"SwitchStatement": function(path){
		var { node } = path
		switchSub.push(node.discriminant.name)
	},
})

//Eval functions in order to get var values
traverse(_testCode, {
	"CallExpression": function(path){
		if(!t.isIdentifier(path.node.callee) && !global[path.node.callee.name] || !global[path.node.callee.name]?.toString()) return;
		let tempCode = global[path.node.callee.name].toString()
		try {
			let ast2 = parser.parse(tempCode)
			traverse(ast2, {
				AssignmentExpression: function(path){
					var { left, right } = path.node
					if(switchSub.includes(left.name)) return;
					try {
						// +-* something that doesnt exists should be prevented
						if(["*=", "+=", "-=", "/="].includes(path.node.operator) && typeof global[left.name] == "undefined") return;
						eval("global." + tempCode.slice(path.node.start, path.node.end))
					}catch(e){
					}
				}
			})
		}catch(e){}
		
	}
})

traverse(_testCode, {
	//Change var xxx; -> var xxx = value;
	"VariableDeclarator": function(path){
		if(typeof global[path.node.id.name] == "number"){
			if(path.node.id.name == "wqx") console.log(global[path.node.id.name])
			path.node.init = t.NumericLiteral(global[path.node.id.name])
		}
	},
	//Change xxx = a + b + c; -> xxx = value
	"AssignmentExpression": function(path){
		//issue here ?
		// return;
		var { left, right } = path.node
		
		// Remove useless xxx = a + b + c
		if(typeof global[left.name] == "number" && !switchSub.includes(left.name)){
			try {
				path.remove()
				return;
			}catch{}
			
		}
		
		// if(global[left.name] && t.isFunctionExpression(right)){
			// console.log(left.name, global[left.name].toString())
			// global[left.name] = 
		// }
		
		if(left.type != "Identifier" || right.type != "Identifier") return;
		if(!globalVars().includes(right.name)) return;
		try {
			// xxx -+= AC -> xxx -+= val
			Object.assign(right, t.numericLiteral(global[right.name]))
		}
		catch{}
	},
	//Change case xxx: -> case 69420:
	"SwitchStatement": function(path){
		var { node } = path
		node.cases = node.cases.map(x => {
			if(!isGlobalVars(x.test?.name)) return x;
			else {
				x.test = t.numericLiteral(global[x.test.name])
				return x
			}
		})
	},
	//Do the main switch loop bullshit idek
	"VariableDeclaration": function(path){
		try {
			global[path.node.id.name] = eval(original.slice(path.node.start, path.node.end))
		}catch{
		}
		if(t.isVariableDeclaration(path.node) && path.node.declarations.length == 1){
			let { id, init } = path.node.declarations[0]
			if(id?.name != callingFunc) return;
			
			let callingParamFunc = init.id.name
			let callingParams = init.params.map(x => x.name)
			
			let whileLoopInMain = init.body.body.find(x => x.type == "WhileStatement")
			let whileCond = original.slice(whileLoopInMain.test.loc.start.index, whileLoopInMain.test.loc.end.index)
			
			function loopFunction(startVal){
				console.log("Starting loop", eval(startVal))
				eval(`var ${callingParams[0]} = ${startVal}`)
				while(eval(whileCond)){
					let cases = whileLoopInMain.body.body[0].cases
					cases = cases.find(x => {
						if(t.isIdentifier(x.test)){
							return eval(x.test?.name) == eval(callingParams[0])
						}else{
							return x.test?.value == eval(callingParams[0])
						}
					})
					let { consequent } = cases
					
					consequent = consequent.find(x => t.isBlockStatement(x))
					
					let { body } = consequent
					
					body = body.filter(x => t.isExpressionStatement(x))
					for(let content of body){
						let { expression } = content	
						//Execute the same loop
						if(t.isCallExpression(expression) && expression.callee.name == callingParamFunc){
							loopFunction(eval(expression.arguments.find(x => t.isIdentifier(x)).name))
							
						}
						else if(t.isAssignmentExpression(expression)){
							if(expression.left.name == callingParams[0]){
								eval(original.slice(expression.start, expression.end))
								console.log(`Loop ${eval(startVal)}: ${eval(callingParams[0])}`)
							}else if(expression.operator == "="){
								try {
									eval("global." + original.slice(expression.start, expression.end))
									if(typeof global[expression.left.name] == "number"){
										expression.right = t.NumericLiteral(global[expression.left.name])
									}
								}
								catch(e){
								}
							}else{
								try {
									eval(original.slice(expression.start, expression.end))
								}
								catch{}
							}
						}
						else if(t.isCallExpression(expression) && expression.callee.object?.name == callingParamFunc){
							let loopVal = expression.arguments.find(x => t.isIdentifier(x))
							if(loopVal){
								loopFunction(eval(expression.arguments.find(x => t.isIdentifier(x)).name))
							}
						}
					}
					if(!eval(whileCond)){
						console.log("Exiting loop", eval(startVal))
					}
				}
			}
			loopFunction(callingVal)
		}
	},
	// funct(xxxx,xxxx,....); -> funct(val, val, .....);
	"CallExpression": function(path){
		for(let arg of path.node.arguments){
			if(!t.isIdentifier(arg)) continue;
			if(typeof global[arg.name] != "number") continue;
			Object.assign(arg, t.NumericLiteral(global[arg.name]))
		}
	},
	"ArrayExpression": function(path){
		for(let elem of path.node.elements){
			if(!t.isIdentifier(elem)) continue;
			if(typeof global[elem.name] != "number") continue;
			Object.assign(elem, t.NumericLiteral(global[elem.name]))
		}
	},
	// arr[xxx] => arr[value]
	"MemberExpression": function(path){
		if(typeof global[path.node.property.name] == "number") path.node.property = (t.NumericLiteral(global[path.node.property.name]))
	}
})

traverse(_testCode, {
	"VariableDeclaration": function(path){
		for(let declaration of path.node.declarations){
			if(!t.isFunctionExpression(declaration.init)) continue;
			try {
				eval(`global.${declaration.id.name} = ${original.slice(declaration.init.start, declaration.init.end)}`)
				if(declaration.init.id){
					// var xx = function yy(){} -> yy = xx
					global[declaration.init.id.name] = eval(declaration.id.name)
				}
				
				if(declaration.init.body.body.length == 1){
					if(declaration.init.body.body[0].expression.type == "AssignmentExpression"){
						eval("global." + original.slice(declaration.init.body.body[0].expression.start, declaration.init.body.body[0].expression.end))
					}
				}
			}
			catch(e){
				// console.log(e)
			}
		}
	},
	// Get all the switch subjects so we dont accidentally remove the assignment
	"SwitchStatement": function(path){
		var { node } = path
		switchSub.push(node.discriminant.name)
	},
})
 
// var fx = global;
// console.log(RB(0, 0))
// console.log(RB(186, 692))

// console.log(global)

var mainFunctions = []
for(let i in global){
	if((typeof global[i] == "function") && (global[i].name != "")){
		if((i != global[i].name)){
			if(i != callingFunc) mainFunctions.push(i)
		}
	}
}

console.log(mainFunctions)

let stringFunction = ""
for(let i in global){
	if((typeof global[i] == "function") && (global[i].name == "")){
		if(i != global[i].name){
			for(let j of mainFunctions){
				if(global[i].toString().includes(j + ".apply")){
					let funcAssignCount = [...original.matchAll(i + "=")].length
					if(funcAssignCount != 1) stringFunction = i
				} 
			}
		}
	}
}

let funcRegex = new RegExp(`${stringFunction}=function.*?}`, "g")
let funcTypes = original.match(funcRegex)
var first = funcTypes.find(x => x.includes(callingFunc))
var second = funcTypes.find(x => !x.includes(callingFunc))

traverse(_testCode, {
	CallExpression: function(path){
		if(!t.isMemberExpression(path.node.callee)) return;
		var { callee } = path.node
		if(callee.object.name != objName) return;
		let argLen = path.node.arguments.length
		let str = original.slice(path.node.start, path.node.end)
		try {
			if(argLen == 4){
				// eval("xx." + callee.property.name + "=" + stringFunction)
				// eval("global." + second)
				// var val = eval(str)
				// console.log(val)
				// if(typeof val == "string"){
					// path.replaceWith(t.StringLiteral(val))
				// }
			}
			else{
				eval("xx." + callee.property.name + "=" + stringFunction)
				eval("global." + first)
				var val = eval(str)
				console.log(val)
				if(typeof val == "string"){
					path.replaceWith(t.StringLiteral(val))
				}
			}
		}catch(e){
			console.log(e)
		}
	},
	// "MemberExpression": function(path){
		// if(typeof global[path.node.property.name] == "number") path.node.property = (t.NumericLiteral(global[path.node.property.name]))
	// }
})

original = generator(_testCode, {
	minified: true
}, original).code


var beautifiedCode = beautify(original, {
	indent_size: "\t",
	space_in_empty_paren: true,
	unescape_strings: true
})

fs.writeFileSync("./deobed.js", original)
fs.writeFileSync("./beautified.js", beautifiedCode)