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


console.log(callingFunc, callingVal)

var _testCode = parser.parse(original)
traverse(_testCode, {
	//Eval all functions make it global
	"FunctionDeclaration": function(path){
		try {
			if(t.isFunctionDeclaration(path)) eval(`global["${path.node.id.name}"] = ${original.slice(path.node.start, path.node.end)}`)
			else{
			}
		}catch(e){
		}
	}
})
traverse(_testCode, {
	"VariableDeclaration": function(path){
		for(let declaration of path.node.declarations){
			if(!t.isFunctionExpression(declaration.init)) continue;
			if(declaration.id.name == callingFunc) continue;
			try {
				eval(`global.${declaration.id.name} = ${original.slice(declaration.init.start, declaration.init.end)}`)
			}
			catch(e){
				// console.log(e)
			}
		}
	}
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
					try {
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
			path.node.init = t.NumericLiteral(global[path.node.id.name])
		}
	},
	//Change xxx = a + b + c; -> xxx = value
	"AssignmentExpression": function(path){
		var { left, right } = path.node
		
		//Remove useless xxx = a + b + c
		if(typeof global[left.name] == "number"){
			path.remove()
			return;
		}
		
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
					}
					if(!eval(whileCond)){
						console.log("Exiting loop", eval(startVal))
					}
				}
			}
			loopFunction(callingVal)
		}
	}
})

// console.log(global)
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